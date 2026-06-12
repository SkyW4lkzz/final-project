#include "lbm.hpp"
#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace lbm {
namespace {

inline int index_T(int x, int y, int i, int nx) {
    return (y * nx + x) * Q_T + i;
}

// --- 新增：將對流場的資料輸出為 CSV ---
void write_convection_field(const std::vector<double>& f, const std::vector<double>& T_dist, 
                            const Config& cfg, const std::string& path_str) {
    std::filesystem::path path(path_str);
    if (path.has_parent_path()) {
        std::filesystem::create_directories(path.parent_path());
    }

    std::ofstream out(path);
    if (!out) throw std::runtime_error("failed to open file: " + path_str);

    out << "x,y,rho,ux,uy,T\n";
    out << std::setprecision(6);
    
    const double T0 = 0.5 * (cfg.Th + cfg.Tc);
    for (int y = 0; y < cfg.ny; ++y) {
        for (int x = 0; x < cfg.nx; ++x) {
            double T_macro = 0.0;
            for (int i = 0; i < Q_T; ++i) {
                T_macro += T_dist[index_T(x, y, i, cfg.nx)];
            }
            const double force_y = cfg.gravity * cfg.beta * (T_macro - T0);
            const CellState state = macroscopic(f, x, y, cfg.nx, 0.0, force_y);
            
            out << x << ',' << y << ',' << state.rho << ',' << state.ux << ',' << state.uy << ',' << T_macro << '\n';
        }
    }
}

void collide_and_stream_convection(
    const std::vector<double>& f, std::vector<double>& f_next,
    const std::vector<double>& T_dist, std::vector<double>& T_next,
    const Config& cfg)
{
    std::fill(f_next.begin(), f_next.end(), 0.0);
    std::fill(T_next.begin(), T_next.end(), 0.0);

    const double T0 = 0.5 * (cfg.Th + cfg.Tc);
    const double tau_T = 2.0 * cfg.kappa + 0.5;

#ifdef _OPENMP
#pragma omp parallel for collapse(2)
#endif
    for (int y = 0; y < cfg.ny; ++y) {
        for (int x = 0; x < cfg.nx; ++x) {
            double T_macro = 0.0;
            for (int i = 0; i < Q_T; ++i) {
                T_macro += T_dist[index_T(x, y, i, cfg.nx)];
            }

            const double force_y = cfg.gravity * cfg.beta * (T_macro - T0);
            const double force_x = 0.0; 
            const CellState state = macroscopic(f, x, y, cfg.nx, force_x, force_y);
            
            // 流場 D2Q9
            for (int i = 0; i < Q; ++i) {
                double post_f = 0.0;
                if (y == 0 || y == cfg.ny - 1) { 
                    post_f = f[index(x, y, opposite[i], cfg.nx)];
                } else {
                    const double feq = equilibrium(i, state.rho, state.ux, state.uy);
                    const double forced = forcing_term(i, state.ux, state.uy, force_x, force_y, relaxation_time(cfg));
                    post_f = f[index(x, y, i, cfg.nx)] - (f[index(x, y, i, cfg.nx)] - feq) / relaxation_time(cfg) + forced;
                }

                int nx_pos = (x + cx[i] + cfg.nx) % cfg.nx; 
                int ny_pos = y + cy[i];
                if (ny_pos >= 0 && ny_pos < cfg.ny) {
                    f_next[index(nx_pos, ny_pos, i, cfg.nx)] = post_f;
                }
            }

            // 溫度場 D2Q4
            for (int i = 0; i < Q_T; ++i) {
                double post_T = 0.0;
                const double Teq = 0.25 * T_macro * (1.0 + 2.0 * (cx_T[i] * state.ux + cy_T[i] * state.uy));
                
                if (y > 0 && y < cfg.ny - 1) {
                    post_T = T_dist[index_T(x, y, i, cfg.nx)] - (T_dist[index_T(x, y, i, cfg.nx)] - Teq) / tau_T;
                } else {
                    post_T = T_dist[index_T(x, y, i, cfg.nx)]; 
                }

                int nx_pos = (x + cx_T[i] + cfg.nx) % cfg.nx;
                int ny_pos = y + cy_T[i];
                if (ny_pos >= 0 && ny_pos < cfg.ny) {
                    T_next[index_T(nx_pos, ny_pos, i, cfg.nx)] = post_T;
                }
            }
        }
    }

#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int x = 0; x < cfg.nx; ++x) {
        T_next[index_T(x, 0, 1, cfg.nx)] = cfg.Th 
            - T_next[index_T(x, 0, 0, cfg.nx)] 
            - T_next[index_T(x, 0, 2, cfg.nx)] 
            - T_next[index_T(x, 0, 3, cfg.nx)];

        T_next[index_T(x, cfg.ny - 1, 3, cfg.nx)] = cfg.Tc 
            - T_next[index_T(x, cfg.ny - 1, 0, cfg.nx)] 
            - T_next[index_T(x, cfg.ny - 1, 1, cfg.nx)] 
            - T_next[index_T(x, cfg.ny - 1, 2, cfg.nx)];
    }
}

} // namespace

void run_convection(const Config& cfg) {
    std::vector<char> solid(static_cast<std::size_t>(cfg.nx) * cfg.ny, 0); 
    std::vector<double> f(static_cast<std::size_t>(cfg.nx) * cfg.ny * Q);
    std::vector<double> f_next(f.size());
    initialize(f, cfg, solid);

    std::vector<double> T_dist(static_cast<std::size_t>(cfg.nx) * cfg.ny * Q_T);
    std::vector<double> T_next(T_dist.size());
    
    srand(42); 
    for (int y = 0; y < cfg.ny; ++y) {
        for (int x = 0; x < cfg.nx; ++x) {
            double base_T = cfg.Th - static_cast<double>(y) / (cfg.ny - 1) * (cfg.Th - cfg.Tc);
            if (y > 0 && y < cfg.ny - 1) {
                double noise = ((std::rand() / (double)RAND_MAX) - 0.5) * 0.01 * (cfg.Th - cfg.Tc);
                base_T += noise;
            }
            for (int i = 0; i < Q_T; ++i) {
                T_dist[index_T(x, y, i, cfg.nx)] = 0.25 * base_T; 
            }
        }
    }

    std::cout << "\nD2Q9+D2Q4 Rayleigh-Bénard Convection\n"
              << "nx= " << cfg.nx << ", ny= " << cfg.ny << ", Steps= " << cfg.steps 
              << "\nnu= " << viscosity(cfg) << ", kappa= " << cfg.kappa << ", gravity= " << cfg.gravity
              << "\nTh= " << cfg.Th << ", Tc= " << cfg.Tc << "\n\n";

    const auto start = std::chrono::steady_clock::now();
    for (int step = 1; step <= cfg.steps; ++step) {
        collide_and_stream_convection(f, f_next, T_dist, T_next, cfg);
        f.swap(f_next);
        T_dist.swap(T_next);

        if (step % cfg.report_interval == 0 || step == cfg.steps) {
            std::cout << "Step=" << std::setw(6) << step << " completed.\n";
        }
        
        // --- 新增：儲存快照 ---
        if (cfg.snapshot_interval > 0 && (step % cfg.snapshot_interval == 0 || step == cfg.steps)) {
            write_convection_field(f, T_dist, cfg, field_snapshot_path(cfg, step));
        }
    }
    const auto stop = std::chrono::steady_clock::now();
    
    // --- 新增：輸出最終結果 ---
    write_convection_field(f, T_dist, cfg, cfg.output);

    const double seconds = std::chrono::duration<double>(stop - start).count();
    std::cout << "\nRuntime= " << seconds << " s\nOutput saved to: " << cfg.output << "\n";
}

} // namespace lbm