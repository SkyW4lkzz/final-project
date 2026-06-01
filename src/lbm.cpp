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
    // Distribution index in a flattened node-major array
    int index(const int x, const int y, const int i, const int nx)
    {
        return ((y * nx + x) * Q) + i;
    }
    // Node index in a flattened 2D array
    int node_index(const int x, const int y, const int nx)
    {
        return y * nx + x;
    }
    // Top and bottom wall test
    bool is_solid(const int y, const int ny)
    {
        return y == 0 || y == ny - 1;
    }
    // Solid-mask lookup
    bool is_solid_node(const std::vector<char> &solid, const int x, const int y, const int nx)
    {
        return solid[node_index(x, y, nx)] != 0;
    }
    // Cylinder center x-position
    int cylinder_x(const Config &cfg)
    {
        return cfg.cylinder_x >= 0 ? cfg.cylinder_x : cfg.nx / 4;
    }
    // Cylinder center y-position
    int cylinder_y(const Config &cfg)
    {
        return cfg.cylinder_y >= 0 ? cfg.cylinder_y : cfg.ny / 2;
    }
    // Cylinder radius from explicit value or grid-size ratio
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
    // Airfoil leading-edge x-position
    int airfoil_x(const Config &cfg)
    {
        return cfg.airfoil_x >= 0 ? cfg.airfoil_x : cfg.nx / 4;
    }
    // Airfoil centerline y-position
    int airfoil_y(const Config &cfg)
    {
        return cfg.airfoil_y >= 0 ? cfg.airfoil_y : cfg.ny / 2;
    }
    // Airfoil chord from explicit value or grid-size ratio
    int airfoil_chord(const Config &cfg)
    {
        if (cfg.airfoil_chord > 0)
        {
            return cfg.airfoil_chord;
        }
        const int reference_size = std::min(cfg.nx, cfg.ny);
        const int chord = static_cast<int>(std::round(cfg.airfoil_chord_ratio * static_cast<double>(reference_size)));
        return std::max(8, chord);
    }
    // Kinematic viscosity
    double viscosity(const Config &cfg)
    {
        return cfg.nu;
    }
    // BGK relaxation time
    double relaxation_time(const Config &cfg)
    {
        return 0.5 + cfg.nu / cs2;
    }
    // Prescribed cylinder inlet velocity
    double cylinder_inlet_ux(const Config &cfg)
    {
        return cfg.inlet_ux;
    }
    // Cylinder viscosity inferred from Reynolds number
    double cylinder_viscosity_from_re(const Config &cfg)
    {
        const double diameter = 2.0 * static_cast<double>(cylinder_radius(cfg));
        return cylinder_inlet_ux(cfg) * diameter / cfg.reynolds_number;
    }
    // Airfoil viscosity inferred from Reynolds number
    double airfoil_viscosity_from_re(const Config &cfg)
    {
        return cylinder_inlet_ux(cfg) * static_cast<double>(airfoil_chord(cfg)) / cfg.reynolds_number;
    }
    // Boundary condition name for output
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
    // Case name for output
    std::string case_name(const CaseType case_type)
    {
        switch (case_type)
        {
        case CaseType::Poiseuille:
            return "poiseuille";
        case CaseType::Cylinder:
            return "cylinder";
        case CaseType::Airfoil:
            return "airfoil";
        }
        return "unknown";
    }
    // Configure OpenMP thread count
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
    // Active OpenMP thread limit
    int max_thread_count()
    {
#ifdef _OPENMP
        return omp_get_max_threads();
#else
        return 1;
#endif
    }
    // D2Q9 equilibrium distribution
    double equilibrium(const int i, const double rho, const double ux, const double uy)
    {
        const double cu = static_cast<double>(cx[i]) * ux + static_cast<double>(cy[i]) * uy;
        const double u2 = ux * ux + uy * uy;
        return w[i] * rho * (1.0 + cu / cs2 + 0.5 * cu * cu / cs4 - 0.5 * u2 / cs2);
    }
    // Guo forcing term for body-force driven flow
    double forcing_term(const int i, const double ux, const double uy, const double force_x, const double force_y, const double tau)
    {
        const double ci_dot_u = static_cast<double>(cx[i]) * ux + static_cast<double>(cy[i]) * uy;
        const double term_x = (static_cast<double>(cx[i]) - ux) / cs2 + static_cast<double>(cx[i]) * ci_dot_u / cs4;
        const double term_y = (static_cast<double>(cy[i]) - uy) / cs2 + static_cast<double>(cy[i]) * ci_dot_u / cs4;
        return w[i] * (1.0 - 0.5 / tau) * (term_x * force_x + term_y * force_y);
    }
    // Macroscopic density and velocity
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
    // Solid mask for walls and optional obstacle
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
        else if (cfg.case_type == CaseType::Airfoil)
        {
            const int x0 = airfoil_x(cfg);
            const int y0 = airfoil_y(cfg);
            const int chord = airfoil_chord(cfg);
            const double thickness = cfg.airfoil_thickness;
            const double angle = -cfg.airfoil_angle * std::acos(-1.0) / 180.0;
            const double cos_angle = std::cos(angle);
            const double sin_angle = std::sin(angle);
            const double rotation_x = static_cast<double>(x0) + 0.25 * static_cast<double>(chord);
            const double rotation_y = static_cast<double>(y0);

#ifdef _OPENMP
#pragma omp parallel for collapse(2)
#endif
            for (int y = 1; y < cfg.ny - 1; ++y)
            {
                for (int x = 1; x < cfg.nx - 1; ++x)
                {
                    const double dx = static_cast<double>(x) - rotation_x;
                    const double dy = static_cast<double>(y) - rotation_y;
                    const double local_x = rotation_x + dx * cos_angle + dy * sin_angle;
                    const double local_y = rotation_y - dx * sin_angle + dy * cos_angle;

                    const double xi = (local_x - static_cast<double>(x0)) / static_cast<double>(chord);
                    if (xi < 0.0 || xi > 1.0)
                    {
                        continue;
                    }

                    const double yt = 5.0 * thickness * static_cast<double>(chord) *
                                      (0.2969 * std::sqrt(xi) - 0.1260 * xi - 0.3516 * xi * xi +
                                       0.2843 * xi * xi * xi - 0.1036 * xi * xi * xi * xi);
                    if (std::abs(local_y - static_cast<double>(y0)) <= std::max(0.5, yt))
                    {
                        solid[node_index(x, y, cfg.nx)] = 1;
                    }
                }
            }
        }

        return solid;
    }
    // Initialize all distributions at equilibrium
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
                const bool inlet_flow = (cfg.case_type == CaseType::Cylinder || cfg.case_type == CaseType::Airfoil) && fluid;
                const double ux = inlet_flow ? cylinder_inlet_ux(cfg) : 0.0;
                const double uy = 0.0;
                for (int i = 0; i < Q; ++i)
                {
                    f[index(x, y, i, cfg.nx)] = equilibrium(i, cfg.outlet_rho, ux, uy);
                }
            }
        }
    }
    // Poiseuille profile snapshot path
    std::string profile_snapshot_path(const Config &cfg, const int step)
    {
        std::ostringstream filename;
        filename << "profile_" << std::setw(7) << std::setfill('0') << step << ".csv";
        return (std::filesystem::path(cfg.snapshot_dir) / filename.str()).string();
    }
    // Velocity field snapshot path
    std::string field_snapshot_path(const Config &cfg, const int step)
    {
        std::ostringstream filename;
        filename << "field_" << std::setw(7) << std::setfill('0') << step << ".csv";
        return (std::filesystem::path(cfg.snapshot_dir) / filename.str()).string();
    }
    // Write a full 2D field CSV
    void write_field(const std::vector<double> &f, const Config &cfg, const std::vector<char> &solid, const std::string &output_path)
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
