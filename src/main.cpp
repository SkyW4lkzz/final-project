#include <array>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace
{
    constexpr int Q = 9;
    constexpr double cs2 = 1.0 / 3.0;
    constexpr double cs4 = cs2 * cs2;
    constexpr std::array<int, Q> cx{0, 1, 0, -1, 0, 1, -1, -1, 1};
    constexpr std::array<int, Q> cy{0, 0, 1, 0, -1, 1, 1, -1, -1};
    constexpr std::array<int, Q> opposite{0, 3, 4, 1, 2, 7, 8, 5, 6};
    constexpr std::array<double, Q> w{
        4.0 / 9.0,
        1.0 / 9.0,
        1.0 / 9.0,
        1.0 / 9.0,
        1.0 / 9.0,
        1.0 / 36.0,
        1.0 / 36.0,
        1.0 / 36.0,
        1.0 / 36.0,
    };

    enum class BoundaryCondition
    {
        OnGrid,
        MidGrid,
    };

    struct Config
    {
        int nx = 40;
        int ny = 33;
        int steps = 10000;
        int report_interval = 1000;
        int snapshot_interval = 0;
        double tau = 1.0;
        double force_x = 0.00125;
        BoundaryCondition boundary = BoundaryCondition::OnGrid;
        std::string output = "test/profile.csv";
        std::string snapshot_dir = "test/snapshots";
    };

    struct CellState
    {
        double rho = 1.0;
        double ux = 0.0;
        double uy = 0.0;
    };

    int index(const int x, const int y, const int i, const int nx)
    {
        return ((y * nx + x) * Q) + i;
    }

    bool is_solid(const int y, const int ny)
    {
        return y == 0 || y == ny - 1;
    }

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

    double equilibrium(const int i, const double rho, const double ux, const double uy)
    {
        const double cu = static_cast<double>(cx[i]) * ux + static_cast<double>(cy[i]) * uy;
        const double u2 = ux * ux + uy * uy;
        return w[i] * rho * (1.0 + cu / cs2 + 0.5 * cu * cu / cs4 - 0.5 * u2 / cs2);
    }

    double forcing_term(const int i, const double ux, const double uy, const double force_x, const double force_y,
                        const double tau)
    {
        const double ci_dot_u = static_cast<double>(cx[i]) * ux + static_cast<double>(cy[i]) * uy;
        const double term_x = (static_cast<double>(cx[i]) - ux) / cs2 + static_cast<double>(cx[i]) * ci_dot_u / cs4;
        const double term_y = (static_cast<double>(cy[i]) - uy) / cs2 + static_cast<double>(cy[i]) * ci_dot_u / cs4;
        return w[i] * (1.0 - 0.5 / tau) * (term_x * force_x + term_y * force_y);
    }

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

    void initialize(std::vector<double> &f, const Config &cfg)
    {
        for (int y = 0; y < cfg.ny; ++y)
        {
            for (int x = 0; x < cfg.nx; ++x)
            {
                for (int i = 0; i < Q; ++i)
                {
                    f[index(x, y, i, cfg.nx)] = equilibrium(i, 1.0, 0.0, 0.0);
                }
            }
        }
    }

    void collide_and_stream(const std::vector<double> &f, std::vector<double> &f_next, const Config &cfg)
    {
        std::fill(f_next.begin(), f_next.end(), 0.0);

        if (cfg.boundary == BoundaryCondition::OnGrid)
        {
            for (int y = 0; y < cfg.ny; ++y)
            {
                for (int x = 0; x < cfg.nx; ++x)
                {
                    for (int i = 0; i < Q; ++i)
                    {
                        double post_collision = 0.0;
                        if (is_solid(y, cfg.ny))
                        {
                            post_collision = f[index(x, y, opposite[i], cfg.nx)];
                        }
                        else
                        {
                            const CellState state = macroscopic(f, x, y, cfg.nx, cfg.force_x, 0.0);
                            const double feq = equilibrium(i, state.rho, state.ux, state.uy);
                            const double forced = forcing_term(i, state.ux, state.uy, cfg.force_x, 0.0, cfg.tau);
                            post_collision =
                                f[index(x, y, i, cfg.nx)] - (f[index(x, y, i, cfg.nx)] - feq) / cfg.tau + forced;
                        }

                        const int nx = (x + cx[i] + cfg.nx) % cfg.nx;
                        const int ny = y + cy[i];
                        if (ny >= 0 && ny < cfg.ny)
                        {
                            f_next[index(nx, ny, i, cfg.nx)] += post_collision;
                        }
                    }
                }
            }
            return;
        }

        for (int y = 1; y < cfg.ny - 1; ++y)
        {
            for (int x = 0; x < cfg.nx; ++x)
            {
                const CellState state = macroscopic(f, x, y, cfg.nx, cfg.force_x, 0.0);

                for (int i = 0; i < Q; ++i)
                {
                    const double feq = equilibrium(i, state.rho, state.ux, state.uy);
                    const double forced = forcing_term(i, state.ux, state.uy, cfg.force_x, 0.0, cfg.tau);
                    const double post_collision =
                        f[index(x, y, i, cfg.nx)] - (f[index(x, y, i, cfg.nx)] - feq) / cfg.tau + forced;

                    const int nx = (x + cx[i] + cfg.nx) % cfg.nx;
                    const int ny = y + cy[i];
                    if (is_solid(ny, cfg.ny))
                    {
                        f_next[index(x, y, opposite[i], cfg.nx)] += post_collision;
                    }
                    else
                    {
                        f_next[index(nx, ny, i, cfg.nx)] += post_collision;
                    }
                }
            }
        }
    }

    double analytical_poiseuille_ux(const int y, const Config &cfg)
    {
        const double nu = cs2 * (cfg.tau - 0.5);
        const double height =
            cfg.boundary == BoundaryCondition::OnGrid ? static_cast<double>(cfg.ny - 1) : static_cast<double>(cfg.ny - 2);
        const double yp = cfg.boundary == BoundaryCondition::OnGrid ? static_cast<double>(y) : static_cast<double>(y) - 0.5;
        return cfg.force_x * yp * (height - yp) / (2.0 * nu);
    }

    struct ErrorMetrics
    {
        double l2_relative = 0.0;
        double max_absolute = 0.0;
    };

    ErrorMetrics compute_error(const std::vector<double> &f, const Config &cfg)
    {
        double numerator = 0.0;
        double denominator = 0.0;
        double max_absolute = 0.0;

        for (int y = 1; y < cfg.ny - 1; ++y)
        {
            double ux_average = 0.0;
            for (int x = 0; x < cfg.nx; ++x)
            {
                ux_average += macroscopic(f, x, y, cfg.nx, cfg.force_x, 0.0).ux;
            }
            ux_average /= static_cast<double>(cfg.nx);

            const double exact = analytical_poiseuille_ux(y, cfg);
            const double diff = ux_average - exact;
            numerator += diff * diff;
            denominator += exact * exact;
            max_absolute = std::max(max_absolute, std::abs(diff));
        }

        const double relative = denominator > 0.0 ? std::sqrt(numerator / denominator) : std::sqrt(numerator);
        return {relative, max_absolute};
    }

    void write_profile(const std::vector<double> &f, const Config &cfg, const std::string &output_path)
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

        out << "y,ux_average,ux_analytical,error\n";
        out << std::setprecision(12);
        for (int y = 1; y < cfg.ny - 1; ++y)
        {
            double ux_average = 0.0;
            for (int x = 0; x < cfg.nx; ++x)
            {
                ux_average += macroscopic(f, x, y, cfg.nx, cfg.force_x, 0.0).ux;
            }
            ux_average /= static_cast<double>(cfg.nx);

            const double exact = analytical_poiseuille_ux(y, cfg);
            out << y << ',' << ux_average << ',' << exact << ',' << ux_average - exact << '\n';
        }
    }

    std::string snapshot_path(const Config &cfg, const int step)
    {
        std::ostringstream filename;
        filename << "profile_" << std::setw(7) << std::setfill('0') << step << ".csv";
        return (std::filesystem::path(cfg.snapshot_dir) / filename.str()).string();
    }

    std::string field_snapshot_path(const Config &cfg, const int step)
    {
        std::ostringstream filename;
        filename << "field_" << std::setw(7) << std::setfill('0') << step << ".csv";
        return (std::filesystem::path(cfg.snapshot_dir) / filename.str()).string();
    }

    void write_field(const std::vector<double> &f, const Config &cfg, const std::string &output_path)
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

        out << "x,y,rho,ux,uy,speed,solid\n";
        out << std::setprecision(12);
        for (int y = 0; y < cfg.ny; ++y)
        {
            for (int x = 0; x < cfg.nx; ++x)
            {
                if (is_solid(y, cfg.ny))
                {
                    out << x << ',' << y << ",0,0,0,0,1\n";
                    continue;
                }

                const CellState state = macroscopic(f, x, y, cfg.nx, cfg.force_x, 0.0);
                const double speed = std::sqrt(state.ux * state.ux + state.uy * state.uy);
                out << x << ',' << y << ',' << state.rho << ',' << state.ux << ',' << state.uy << ',' << speed
                    << ",0\n";
            }
        }
    }

    int parse_int(const char *value, const std::string &name)
    {
        char *end = nullptr;
        const long parsed = std::strtol(value, &end, 10);
        if (*end != '\0')
        {
            throw std::invalid_argument("invalid integer for " + name + ": " + value);
        }
        return static_cast<int>(parsed);
    }

    double parse_double(const char *value, const std::string &name)
    {
        char *end = nullptr;
        const double parsed = std::strtod(value, &end);
        if (*end != '\0')
        {
            throw std::invalid_argument("invalid floating-point value for " + name + ": " + value);
        }
        return parsed;
    }

    void parse_grid(const std::string &value, Config &cfg)
    {
        const std::size_t separator = value.find_first_of("xX,");
        if (separator == std::string::npos)
        {
            const int size = parse_int(value.c_str(), "--grid");
            cfg.nx = size;
            cfg.ny = size;
            return;
        }

        cfg.nx = parse_int(value.substr(0, separator).c_str(), "--grid nx");
        cfg.ny = parse_int(value.substr(separator + 1).c_str(), "--grid ny");
    }

    BoundaryCondition parse_boundary(const std::string &value)
    {
        if (value == "on-grid" || value == "ongrid")
        {
            return BoundaryCondition::OnGrid;
        }
        if (value == "mid-grid" || value == "midgrid" || value == "half-way" || value == "halfway")
        {
            return BoundaryCondition::MidGrid;
        }
        throw std::invalid_argument("invalid boundary condition: " + value + " (expected on-grid or mid-grid)");
    }

    Config parse_args(const int argc, char **argv)
    {
        Config cfg;
        for (int i = 1; i < argc; ++i)
        {
            const std::string arg = argv[i];
            auto require_value = [&](const std::string &option)
            {
                if (i + 1 >= argc)
                {
                    throw std::invalid_argument("missing value for " + option);
                }
                return argv[++i];
            };

            if (arg == "--grid")
            {
                parse_grid(require_value(arg), cfg);
            }
            else if (arg == "--nx")
            {
                cfg.nx = parse_int(require_value(arg), arg);
            }
            else if (arg == "--ny")
            {
                cfg.ny = parse_int(require_value(arg), arg);
            }
            else if (arg == "--steps")
            {
                cfg.steps = parse_int(require_value(arg), arg);
            }
            else if (arg == "--tau")
            {
                cfg.tau = parse_double(require_value(arg), arg);
            }
            else if (arg == "--force-x")
            {
                cfg.force_x = parse_double(require_value(arg), arg);
            }
            else if (arg == "--boundary")
            {
                cfg.boundary = parse_boundary(require_value(arg));
            }
            else if (arg == "--report-interval")
            {
                cfg.report_interval = parse_int(require_value(arg), arg);
            }
            else if (arg == "--output")
            {
                cfg.output = require_value(arg);
            }
            else if (arg == "--snapshot-interval")
            {
                cfg.snapshot_interval = parse_int(require_value(arg), arg);
            }
            else if (arg == "--snapshot-dir")
            {
                cfg.snapshot_dir = require_value(arg);
            }
            else if (arg == "--help")
            {
                std::cout << "Usage: lbm_solver [--grid N|NxM] [--nx N] [--ny N] [--steps N] [--tau T]\n"
                          << "                  [--force-x F] [--boundary on-grid|mid-grid]\n"
                          << "                  [--report-interval N] [--output path]\n"
                          << "                  [--snapshot-interval N] [--snapshot-dir path]\n";
                std::exit(0);
            }
            else
            {
                throw std::invalid_argument("unknown option: " + arg);
            }
        }

        if (cfg.nx < 3 || cfg.ny < 4)
        {
            throw std::invalid_argument("Domain must satisfy nx >= 3 and ny >= 4");
        }
        if (cfg.steps < 0 || cfg.report_interval < 1)
        {
            throw std::invalid_argument("Steps must be non-negative and report interval must be positive");
        }
        if (cfg.snapshot_interval < 0)
        {
            throw std::invalid_argument("snapshot interval must be non-negative");
        }
        if (cfg.tau <= 0.5)
        {
            throw std::invalid_argument("tau must be greater than 0.5 for positive viscosity");
        }
        return cfg;
    }

} // namespace

int main(int argc, char **argv)
{
    try
    {
        const Config cfg = parse_args(argc, argv);
        std::vector<double> f(static_cast<std::size_t>(cfg.nx) * cfg.ny * Q);
        std::vector<double> f_next(f.size());
        initialize(f, cfg);

        std::cout << "D2Q9 BGK Poiseuille flow\n"
                  << "nx=" << cfg.nx << " ny=" << cfg.ny << " steps=" << cfg.steps
                  << " tau=" << cfg.tau << " force_x=" << cfg.force_x
                  << " boundary=" << boundary_name(cfg.boundary) << '\n';

        if (cfg.snapshot_interval > 0)
        {
            write_profile(f, cfg, snapshot_path(cfg, 0));
            write_field(f, cfg, field_snapshot_path(cfg, 0));
        }

        const auto start = std::chrono::steady_clock::now();
        for (int step = 1; step <= cfg.steps; ++step)
        {
            collide_and_stream(f, f_next, cfg);
            f.swap(f_next);

            if (step % cfg.report_interval == 0 || step == cfg.steps)
            {
                const ErrorMetrics error = compute_error(f, cfg);
                std::cout << "step=" << std::setw(7) << step << " relative_l2=" << std::scientific
                          << error.l2_relative << " max_abs=" << error.max_absolute << std::defaultfloat
                          << '\n';
            }

            if (cfg.snapshot_interval > 0 && (step % cfg.snapshot_interval == 0 || step == cfg.steps))
            {
                write_profile(f, cfg, snapshot_path(cfg, step));
                write_field(f, cfg, field_snapshot_path(cfg, step));
            }
        }
        const auto stop = std::chrono::steady_clock::now();

        const double seconds = std::chrono::duration<double>(stop - start).count();
        const double fluid_nodes = static_cast<double>(cfg.nx) * static_cast<double>(cfg.ny - 2);
        const double mlups = seconds > 0.0 ? fluid_nodes * static_cast<double>(cfg.steps) / seconds / 1.0e6 : 0.0;
        const ErrorMetrics error = compute_error(f, cfg);

        write_profile(f, cfg, cfg.output);
        std::cout << "finished in " << seconds << " s, " << mlups << " MLUPS\n"
                  << "final relative_l2=" << std::scientific << error.l2_relative
                  << " max_abs=" << error.max_absolute << std::defaultfloat << '\n'
                  << "wrote " << cfg.output << '\n';
    }
    catch (const std::exception &error)
    {
        std::cerr << "error: " << error.what() << '\n';
        return 1;
    }
    return 0;
}
