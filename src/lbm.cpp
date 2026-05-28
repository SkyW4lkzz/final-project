#include "lbm.hpp"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace lbm
{
    //
    int index(const int x, const int y, const int i, const int nx)
    {
        return ((y * nx + x) * Q) + i;
    }
    //
    int node_index(const int x, const int y, const int nx)
    {
        return y * nx + x;
    }
    //
    bool is_solid(const int y, const int ny)
    {
        return y == 0 || y == ny - 1;
    }
    //
    bool is_solid_node(const std::vector<char> &solid, const int x, const int y, const int nx)
    {
        return solid[node_index(x, y, nx)] != 0;
    }
    //
    int cylinder_x(const Config &cfg)
    {
        return cfg.cylinder_x >= 0 ? cfg.cylinder_x : cfg.nx / 4;
    }
    //
    int cylinder_y(const Config &cfg)
    {
        return cfg.cylinder_y >= 0 ? cfg.cylinder_y : cfg.ny / 2;
    }
    //
    int cylinder_radius(const Config &cfg)
    {
        if (cfg.cylinder_radius > 0)
        {
            return cfg.cylinder_radius;
        }
        const int reference_size = std::min(cfg.nx, cfg.ny);
        const int radius = static_cast<int>(std::round(cfg.cylinder_radius_ratio * static_cast<double>(reference_size)));
        return std::max(2, radius);
    }
    //
    double viscosity(const Config &cfg)
    {
        return cfg.nu;
    }
    //
    double relaxation_time(const Config &cfg)
    {
        return 0.5 + cfg.nu / cs2;
    }
    //
    double cylinder_inlet_ux(const Config &cfg)
    {
        const double diameter = 2.0 * static_cast<double>(cylinder_radius(cfg));
        return cfg.reynolds_number * viscosity(cfg) / diameter;
    }
    //
    std::string boundary_name(const BoundaryCondition boundary)
    {
        switch (boundary)
        {
        case BoundaryCondition::OnGrid:
            return "on-grid";
        case BoundaryCondition::MidGrid:
            return "mid-grid";
        }
        return "unknown";
    }
    //
    std::string case_name(const CaseType case_type)
    {
        switch (case_type)
        {
        case CaseType::Poiseuille:
            return "poiseuille";
        case CaseType::Cylinder:
            return "cylinder";
        }
        return "unknown";
    }
    //
    void configure_threads(const Config &cfg)
    {
#ifdef _OPENMP
        if (cfg.threads > 0)
        {
            omp_set_num_threads(cfg.threads);
        }
#else
        if (cfg.threads > 1)
        {
            throw std::invalid_argument("OpenMP is not enabled in this build; use --threads 1 or rebuild with OpenMP");
        }
#endif
    }
    //
    int max_thread_count()
    {
#ifdef _OPENMP
        return omp_get_max_threads();
#else
        return 1;
#endif
    }
    //
    double equilibrium(const int i, const double rho, const double ux, const double uy)
    {
        const double cu = static_cast<double>(cx[i]) * ux + static_cast<double>(cy[i]) * uy;
        const double u2 = ux * ux + uy * uy;
        return w[i] * rho * (1.0 + cu / cs2 + 0.5 * cu * cu / cs4 - 0.5 * u2 / cs2);
    }
    //
    double forcing_term(const int i, const double ux, const double uy, const double force_x, const double force_y,
                        const double tau)
    {
        const double ci_dot_u = static_cast<double>(cx[i]) * ux + static_cast<double>(cy[i]) * uy;
        const double term_x = (static_cast<double>(cx[i]) - ux) / cs2 + static_cast<double>(cx[i]) * ci_dot_u / cs4;
        const double term_y = (static_cast<double>(cy[i]) - uy) / cs2 + static_cast<double>(cy[i]) * ci_dot_u / cs4;
        return w[i] * (1.0 - 0.5 / tau) * (term_x * force_x + term_y * force_y);
    }
    //
    CellState macroscopic(const std::vector<double> &f, const int x, const int y, const int nx, const double force_x,
                          const double force_y)
    {
        CellState state;
        state.rho = 0.0;
        double momentum_x = 0.0;
        double momentum_y = 0.0;

        for (int i = 0; i < Q; ++i)
        {
            const double fi = f[index(x, y, i, nx)];
            state.rho += fi;
            momentum_x += static_cast<double>(cx[i]) * fi;
            momentum_y += static_cast<double>(cy[i]) * fi;
        }

        state.ux = (momentum_x + 0.5 * force_x) / state.rho;
        state.uy = (momentum_y + 0.5 * force_y) / state.rho;
        return state;
    }
    //
    std::vector<char> build_solid_mask(const Config &cfg)
    {
        std::vector<char> solid(static_cast<std::size_t>(cfg.nx) * cfg.ny, 0);

        for (int x = 0; x < cfg.nx; ++x)
        {
            solid[node_index(x, 0, cfg.nx)] = 1;
            solid[node_index(x, cfg.ny - 1, cfg.nx)] = 1;
        }

        if (cfg.case_type == CaseType::Cylinder)
        {
            const int cx0 = cylinder_x(cfg);
            const int cy0 = cylinder_y(cfg);
            const int radius = cylinder_radius(cfg);
            const int radius2 = radius * radius;

#ifdef _OPENMP
#pragma omp parallel for collapse(2)
#endif
            for (int y = 1; y < cfg.ny - 1; ++y)
            {
                for (int x = 1; x < cfg.nx - 1; ++x)
                {
                    const int dx = x - cx0;
                    const int dy = y - cy0;
                    if (dx * dx + dy * dy <= radius2)
                    {
                        solid[node_index(x, y, cfg.nx)] = 1;
                    }
                }
            }
        }

        return solid;
    }
    //
    void initialize(std::vector<double> &f, const Config &cfg, const std::vector<char> &solid)
    {
#ifdef _OPENMP
#pragma omp parallel for collapse(2)
#endif
        for (int y = 0; y < cfg.ny; ++y)
        {
            for (int x = 0; x < cfg.nx; ++x)
            {
                const bool fluid = !is_solid_node(solid, x, y, cfg.nx);
                const double ux = cfg.case_type == CaseType::Cylinder && fluid ? cylinder_inlet_ux(cfg) : 0.0;
                const double uy = 0.0;
                for (int i = 0; i < Q; ++i)
                {
                    f[index(x, y, i, cfg.nx)] = equilibrium(i, cfg.outlet_rho, ux, uy);
                }
            }
        }
    }
    //
    std::string profile_snapshot_path(const Config &cfg, const int step)
    {
        std::ostringstream filename;
        filename << "profile_" << std::setw(7) << std::setfill('0') << step << ".csv";
        return (std::filesystem::path(cfg.snapshot_dir) / filename.str()).string();
    }
    //
    std::string field_snapshot_path(const Config &cfg, const int step)
    {
        std::ostringstream filename;
        filename << "field_" << std::setw(7) << std::setfill('0') << step << ".csv";
        return (std::filesystem::path(cfg.snapshot_dir) / filename.str()).string();
    }
    //
    void write_field(const std::vector<double> &f, const Config &cfg, const std::vector<char> &solid,
                     const std::string &output_path)
    {
        const std::filesystem::path path(output_path);
        if (path.has_parent_path())
        {
            std::filesystem::create_directories(path.parent_path());
        }

        std::ofstream out(path);
        if (!out)
        {
            throw std::runtime_error("failed to open output file: " + output_path);
        }

        out << "x,y,rho,ux,uy,velocity,solid\n";
        out << std::setprecision(12);
        for (int y = 0; y < cfg.ny; ++y)
        {
            for (int x = 0; x < cfg.nx; ++x)
            {
                if (is_solid_node(solid, x, y, cfg.nx))
                {
                    out << x << ',' << y << ",0,0,0,0,1\n";
                    continue;
                }

                const double force_x = cfg.case_type == CaseType::Poiseuille ? cfg.force_x : 0.0;
                const CellState state = macroscopic(f, x, y, cfg.nx, force_x, 0.0);
                const double velocity = std::sqrt(state.ux * state.ux + state.uy * state.uy);
                out << x << ',' << y << ',' << state.rho << ',' << state.ux << ',' << state.uy << ',' << velocity << ",0\n";
            }
        }
    }
} // namespace lbm
