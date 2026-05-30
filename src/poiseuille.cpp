#include "lbm.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace lbm
{
    namespace
    {
        //
        void collide_and_stream_poiseuille(const std::vector<double> &f, std::vector<double> &f_next, const Config &cfg)
        {
            std::fill(f_next.begin(), f_next.end(), 0.0);

            if (cfg.boundary == BoundaryCondition::OnGrid)
            {
#ifdef _OPENMP
#pragma omp parallel for collapse(2)
#endif
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
                                const double forced = forcing_term(i, state.ux, state.uy, cfg.force_x, 0.0, relaxation_time(cfg));
                                post_collision =
                                    f[index(x, y, i, cfg.nx)] - (f[index(x, y, i, cfg.nx)] - feq) / relaxation_time(cfg) + forced;
                            }

                            const int nx = (x + cx[i] + cfg.nx) % cfg.nx;
                            const int ny = y + cy[i];
                            if (ny >= 0 && ny < cfg.ny)
                            {
                                f_next[index(nx, ny, i, cfg.nx)] = post_collision;
                            }
                        }
                    }
                }
                return;
            }

#ifdef _OPENMP
#pragma omp parallel for collapse(2)
#endif
            for (int y = 1; y < cfg.ny - 1; ++y)
            {
                for (int x = 0; x < cfg.nx; ++x)
                {
                    const CellState state = macroscopic(f, x, y, cfg.nx, cfg.force_x, 0.0);

                    for (int i = 0; i < Q; ++i)
                    {
                        const double feq = equilibrium(i, state.rho, state.ux, state.uy);
                        const double forced = forcing_term(i, state.ux, state.uy, cfg.force_x, 0.0, relaxation_time(cfg));
                        const double post_collision =
                            f[index(x, y, i, cfg.nx)] - (f[index(x, y, i, cfg.nx)] - feq) / relaxation_time(cfg) + forced;

                        const int nx = (x + cx[i] + cfg.nx) % cfg.nx;
                        const int ny = y + cy[i];
                        if (is_solid(ny, cfg.ny))
                        {
                            f_next[index(x, y, opposite[i], cfg.nx)] = post_collision;
                        }
                        else
                        {
                            f_next[index(nx, ny, i, cfg.nx)] = post_collision;
                        }
                    }
                }
            }
        }
        //
        double analytical_poiseuille_ux(const int y, const Config &cfg)
        {
            const double nu = viscosity(cfg);
            const double height =
                cfg.boundary == BoundaryCondition::OnGrid ? static_cast<double>(cfg.ny - 1) : static_cast<double>(cfg.ny - 2);
            const double yp = cfg.boundary == BoundaryCondition::OnGrid ? static_cast<double>(y) : static_cast<double>(y) - 0.5;
            return cfg.force_x * yp * (height - yp) / (2.0 * nu);
        }
        //
        ErrorMetrics compute_error(const std::vector<double> &f, const Config &cfg)
        {
            double numerator = 0.0;
            double denominator = 0.0;
            double max_absolute = 0.0;

#ifdef _OPENMP
#pragma omp parallel for reduction(+ : numerator, denominator) reduction(max : max_absolute)
#endif
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
        //
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
            const double height =
                cfg.boundary == BoundaryCondition::OnGrid ? static_cast<double>(cfg.ny - 1) : static_cast<double>(cfg.ny - 2);
            out << 0.0 << ",0,0,0\n";
            for (int y = 1; y < cfg.ny - 1; ++y)
            {
                double ux_average = 0.0;
                for (int x = 0; x < cfg.nx; ++x)
                {
                    ux_average += macroscopic(f, x, y, cfg.nx, cfg.force_x, 0.0).ux;
                }
                ux_average /= static_cast<double>(cfg.nx);

                const double y_phys = cfg.boundary == BoundaryCondition::OnGrid ? static_cast<double>(y) : static_cast<double>(y) - 0.5;
                const double exact = analytical_poiseuille_ux(y, cfg);
                out << y_phys << ',' << ux_average << ',' << exact << ',' << ux_average - exact << '\n';
            }
            out << height << ",0,0,0\n";
        }
        //
        std::string summary_path_for(const std::string &profile_output_path)
        {
            const std::filesystem::path profile_path(profile_output_path);
            const std::filesystem::path parent = profile_path.has_parent_path() ? profile_path.parent_path() : std::filesystem::path(".");
            return (parent / "poiseuille_summary.csv").string();
        }
        //
        void write_summary(const Config &cfg, const double seconds, const double mlups, const ErrorMetrics &error,
                           const std::string &summary_path)
        {
            const std::filesystem::path path(summary_path);
            if (path.has_parent_path())
            {
                std::filesystem::create_directories(path.parent_path());
            }

            const bool write_header = !std::filesystem::exists(path) || std::filesystem::file_size(path) == 0;
            std::ofstream out(path, std::ios::app);
            if (!out)
            {
                throw std::runtime_error("failed to open summary file: " + summary_path);
            }

            if (write_header)
            {
                out << "case,nx,ny,nu,tau,force_x,boundary,steps,threads,computing_time,mlups,relative_l2,profile_output\n";
            }

            out << std::setprecision(12)
                << "poiseuille," << cfg.nx << ',' << cfg.ny << ',' << viscosity(cfg) << ',' << relaxation_time(cfg) << ','
                << cfg.force_x << ',' << boundary_name(cfg.boundary) << ',' << cfg.steps << ',' << max_thread_count() << ','
                << seconds << ',' << mlups << ',' << error.l2_relative << ',' << cfg.output << '\n';
        }
    } // namespace
    //
    void run_poiseuille(const Config &cfg)
    {
        const std::vector<char> solid = build_solid_mask(cfg);
        std::vector<double> f(static_cast<std::size_t>(cfg.nx) * cfg.ny * Q);
        std::vector<double> f_next(f.size());
        initialize(f, cfg, solid);

        std::cout << "\nD2Q9 BGK Plane Poiseuille Flow\n"
                  << "nx= " << cfg.nx << ", ny= " << cfg.ny << ", Total steps= " << cfg.steps << ", nu= " << viscosity(cfg)
                  << ", tau= " << relaxation_time(cfg) << ", Threads= " << max_thread_count() << ", force_x= " << cfg.force_x
                  << ", Boundary= " << boundary_name(cfg.boundary) << "\n\n";

        if (cfg.snapshot_interval > 0)
        {
            write_profile(f, cfg, profile_snapshot_path(cfg, 0));
            write_field(f, cfg, solid, field_snapshot_path(cfg, 0));
        }

        const auto start = std::chrono::steady_clock::now();
        for (int step = 1; step <= cfg.steps; ++step)
        {
            collide_and_stream_poiseuille(f, f_next, cfg);
            f.swap(f_next);

            if (step % cfg.report_interval == 0 || step == cfg.steps)
            {
                const ErrorMetrics error = compute_error(f, cfg);
                std::cout << "Step=" << std::setw(6) << step << ", Relative L2= " << std::scientific
                          << error.l2_relative << std::defaultfloat << "\n";
            }

            if (cfg.snapshot_interval > 0 && (step % cfg.snapshot_interval == 0 || step == cfg.steps))
            {
                write_profile(f, cfg, profile_snapshot_path(cfg, step));
                write_field(f, cfg, solid, field_snapshot_path(cfg, step));
            }
        }
        const auto stop = std::chrono::steady_clock::now();

        const double seconds = std::chrono::duration<double>(stop - start).count();
        const double fluid_nodes = static_cast<double>(std::count(solid.begin(), solid.end(), 0));
        const double mlups = seconds > 0.0 ? fluid_nodes * static_cast<double>(cfg.steps) / seconds / 1.0e6 : 0.0;
        const ErrorMetrics error = compute_error(f, cfg);

        write_profile(f, cfg, cfg.output);
        const std::string summary_path = summary_path_for(cfg.output);
        write_summary(cfg, seconds, mlups, error, summary_path);
        std::cout << "\nRuntime= " << seconds << " s, MLUPS= " << mlups << ", Final Relative L2= " << std::scientific << error.l2_relative << std::defaultfloat << "\n"
                  << "Output= " << cfg.output << "\n"
                  << "Summary= " << summary_path << "\n";
    }
} // namespace lbm
