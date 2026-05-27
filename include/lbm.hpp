#pragma once

#include <array>
#include <string>
#include <vector>

namespace lbm
{
    // D2Q9 constants
    inline constexpr int Q = 9;
    inline constexpr double cs2 = 1.0 / 3.0;
    inline constexpr double cs4 = cs2 * cs2;
    inline constexpr std::array<int, Q> cx{0, 1, 0, -1, 0, 1, -1, -1, 1};
    inline constexpr std::array<int, Q> cy{0, 0, 1, 0, -1, 1, 1, -1, -1};
    inline constexpr std::array<int, Q> opposite{0, 3, 4, 1, 2, 7, 8, 5, 6};
    inline constexpr std::array<double, Q> w{
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

    // Boundary condition
    enum class BoundaryCondition
    {
        OnGrid,
        MidGrid,
    };

    // Plane Poiseuille flow or flow past a cylinder
    enum class CaseType
    {
        Poiseuille,
        Cylinder,
    };

    // Config stores all runtime settings
    struct Config
    {
        int nx = 32;
        int ny = 33;
        int steps = 10000;
        int report_interval = 1000;
        int snapshot_interval = 0;
        int threads = 0;

        double tau = 1.0;
        double force_x = 0.0125 / 32.0;
        double reynolds_number = 40.0;
        double outlet_rho = 1.0;

        int cylinder_x = -1;
        int cylinder_y = -1;
        int cylinder_radius = -1;

        BoundaryCondition boundary = BoundaryCondition::OnGrid;
        CaseType case_type = CaseType::Poiseuille;

        std::string output = "test/poiseuille_profile.csv";
        std::string snapshot_dir = "test/snapshots";
    };

    // CellState stores macroscopic fluid variables at one node
    struct CellState
    {
        double rho = 1.0;
        double ux = 0.0;
        double uy = 0.0;
    };

    // Plane Poiseuille flow use ErrorMetrics to report numerical error against the analytical solution
    struct ErrorMetrics
    {
        double l2_relative = 0.0;
        double max_absolute = 0.0;
    };

    // Utility functions
    int index(int x, int y, int i, int nx);
    int node_index(int x, int y, int nx);

    // Geometry and names
    bool is_solid(int y, int ny);
    bool is_solid_node(const std::vector<char> &solid, int x, int y, int nx);
    int cylinder_x(const Config &cfg);
    int cylinder_y(const Config &cfg);
    int cylinder_radius(const Config &cfg);
    double viscosity(const Config &cfg);
    double cylinder_inlet_ux(const Config &cfg);
    std::string boundary_name(BoundaryCondition boundary);
    std::string case_name(CaseType case_type);

    // OpenMP helpers
    void configure_threads(const Config &cfg);
    int max_thread_count();

    // LBM core functions
    double equilibrium(int i, double rho, double ux, double uy);
    double forcing_term(int i, double ux, double uy, double force_x, double force_y, double tau);
    CellState macroscopic(const std::vector<double> &f, int x, int y, int nx, double force_x, double force_y);

    // Setup and output
    std::vector<char> build_solid_mask(const Config &cfg);
    void initialize(std::vector<double> &f, const Config &cfg, const std::vector<char> &solid);

    //
    std::string profile_snapshot_path(const Config &cfg, int step);
    std::string field_snapshot_path(const Config &cfg, int step);
    void write_field(const std::vector<double> &f, const Config &cfg, const std::vector<char> &solid, const std::string &output_path);

    // Case runners
    void run_poiseuille(const Config &cfg);
    void run_cylinder(const Config &cfg);
}   // namespace lbm
