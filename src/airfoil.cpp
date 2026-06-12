#include "lbm.hpp"

#include <algorithm>
#include <chrono>
#include <iomanip>
#include <iostream>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace lbm
{
    namespace
    {
        // Zou-He velocity boundary at the inlet
        void apply_zou_he_velocity_inlet(std::vector<double> &f, const Config &cfg)
        {
            const int x = 0;
            const double ux = cylinder_inlet_ux(cfg);
#ifdef _OPENMP
#pragma omp parallel for
#endif
            for (int y = 1; y < cfg.ny - 1; ++y)
            {
                const double f0 = f[index(x, y, 0, cfg.nx)];
                const double f2 = f[index(x, y, 2, cfg.nx)];
                const double f3 = f[index(x, y, 3, cfg.nx)];
                const double f4 = f[index(x, y, 4, cfg.nx)];
                const double f6 = f[index(x, y, 6, cfg.nx)];
                const double f7 = f[index(x, y, 7, cfg.nx)];
                const double rho = (f0 + f2 + f4 + 2.0 * (f3 + f6 + f7)) / (1.0 - ux);

                f[index(x, y, 1, cfg.nx)] = f3 + (2.0 / 3.0) * rho * ux;
                f[index(x, y, 5, cfg.nx)] = f7 + 0.5 * (f4 - f2) + (1.0 / 6.0) * rho * ux;
                f[index(x, y, 8, cfg.nx)] = f6 + 0.5 * (f2 - f4) + (1.0 / 6.0) * rho * ux;
            }
        }

        // Zou-He density boundary at the outlet
        void apply_zou_he_pressure_outlet(std::vector<double> &f, const Config &cfg)
        {
            const int x = cfg.nx - 1;
            const double rho = cfg.outlet_rho;
#ifdef _OPENMP
#pragma omp parallel for
#endif
            for (int y = 1; y < cfg.ny - 1; ++y)
            {
                const double f0 = f[index(x, y, 0, cfg.nx)];
                const double f1 = f[index(x, y, 1, cfg.nx)];
                const double f2 = f[index(x, y, 2, cfg.nx)];
                const double f4 = f[index(x, y, 4, cfg.nx)];
                const double f5 = f[index(x, y, 5, cfg.nx)];
                const double f8 = f[index(x, y, 8, cfg.nx)];
                const double ux = -1.0 + (f0 + f2 + f4 + 2.0 * (f1 + f5 + f8)) / rho;

                f[index(x, y, 3, cfg.nx)] = f1 - (2.0 / 3.0) * rho * ux;
                f[index(x, y, 7, cfg.nx)] = f5 + 0.5 * (f2 - f4) - (1.0 / 6.0) * rho * ux;
                f[index(x, y, 6, cfg.nx)] = f8 + 0.5 * (f4 - f2) - (1.0 / 6.0) * rho * ux;
            }
        }

        // Collision, streaming, and bounce-back for airfoil flow
        void collide_and_stream_airfoil(const std::vector<double> &f, std::vector<double> &f_next, const Config &cfg,
                                        const std::vector<char> &solid)
        {
            std::fill(f_next.begin(), f_next.end(), 0.0);

#ifdef _OPENMP
#pragma omp parallel for collapse(2)
#endif
            for (int y = 0; y < cfg.ny; ++y)
            {
                for (int x = 0; x < cfg.nx; ++x)
                {
                    if (is_solid_node(solid, x, y, cfg.nx))
                    {
                        continue;
                    }

                    const CellState state = macroscopic(f, x, y, cfg.nx, 0.0, 0.0);
                    for (int i = 0; i < Q; ++i)
                    {
                        const double feq = equilibrium(i, state.rho, state.ux, state.uy);
                        const double post_collision =
                            f[index(x, y, i, cfg.nx)] - (f[index(x, y, i, cfg.nx)] - feq) / relaxation_time(cfg);

                        const int xn = x + cx[i];
                        const int yn = y + cy[i];
                        if (xn < 0 || xn >= cfg.nx || yn < 0 || yn >= cfg.ny)
                        {
                            continue;
                        }

                        if (is_solid_node(solid, xn, yn, cfg.nx))
                        {
                            f_next[index(x, y, opposite[i], cfg.nx)] = post_collision;
                        }
                        else
                        {
                            f_next[index(xn, yn, i, cfg.nx)] = post_collision;
                        }
                    }
                }
            }

            apply_zou_he_velocity_inlet(f_next, cfg);
            apply_zou_he_pressure_outlet(f_next, cfg);
        }
    } // namespace

    // Run the airfoil case
    void run_airfoil(const Config &cfg)
    {
        const std::vector<char> solid = build_solid_mask(cfg);
        std::vector<double> f(static_cast<std::size_t>(cfg.nx) * cfg.ny * Q);
        std::vector<double> f_next(f.size());
        initialize(f, cfg, solid);

        std::cout << "\nD2Q9 BGK Flow past an Airfoil\n"
                  << "nx= " << cfg.nx << ", ny= " << cfg.ny << ", Total steps= " << cfg.steps << ", nu= " << viscosity(cfg)
                  << ", tau= " << relaxation_time(cfg) << ", Threads= " << max_thread_count() << ", Re= " << cfg.reynolds_number
                  << ", inlet_ux= " << cylinder_inlet_ux(cfg) << ", outlet_rho= " << cfg.outlet_rho
                  << ", Airfoil=(" << airfoil_x(cfg) << "," << airfoil_y(cfg) << "), Chord= " << airfoil_chord(cfg)
                  << ", Thickness= " << cfg.airfoil_thickness << ", Angle= " << cfg.airfoil_angle << " deg\n"
                  << "\nStability range: " << min_cylinder_inlet_ux << " <= inlet_ux <= " << max_cylinder_inlet_ux
                  << ", " << min_cylinder_tau << " <= tau <= " << max_cylinder_tau << "\n\n";

        if (cfg.snapshot_interval > 0)
        {
            write_field(f, cfg, solid, field_snapshot_path(cfg, 0));
        }

        const auto start = std::chrono::steady_clock::now();
        // One collide-stream loop is one LBM time step with dt = 1
        for (int step = 1; step <= cfg.steps; ++step)
        {
            collide_and_stream_airfoil(f, f_next, cfg, solid);
            f.swap(f_next);

            if (step % cfg.report_interval == 0 || step == cfg.steps)
            {
                std::cout << "Step=" << std::setw(6) << step << "\n";
            }

            if (cfg.snapshot_interval > 0 && (step % cfg.snapshot_interval == 0 || step == cfg.steps))
            {
                write_field(f, cfg, solid, field_snapshot_path(cfg, step));
            }
        }
        const auto stop = std::chrono::steady_clock::now();

        const double seconds = std::chrono::duration<double>(stop - start).count();
        const double fluid_nodes = static_cast<double>(std::count(solid.begin(), solid.end(), 0));
        const double mlups = seconds > 0.0 ? fluid_nodes * static_cast<double>(cfg.steps) / seconds / 1.0e6 : 0.0;

        write_field(f, cfg, solid, cfg.output);
        std::cout << "\nRuntime= " << seconds << " s, MLUPS= " << mlups << "\n"
                  << "Output= " << cfg.output << "\n";
    }
} // namespace lbm
