#include "lbm.hpp"

#include <cstdlib>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

namespace
{
    // Parse an integer command-line value
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
    // Parse a floating-point command-line value
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
    // Parse --grid as N or NxM
    bool parse_grid(const std::string &value, lbm::Config &cfg)
    {
        const std::size_t separator = value.find_first_of("xX,");
        if (separator == std::string::npos)
        {
            const int size = parse_int(value.c_str(), "--grid");
            cfg.nx = size;
            cfg.ny = size;
            return true;
        }

        cfg.nx = parse_int(value.substr(0, separator).c_str(), "--grid nx");
        cfg.ny = parse_int(value.substr(separator + 1).c_str(), "--grid ny");
        return false;
    }
    // Convert Poiseuille channel height into lattice rows
    int poiseuille_y_nodes(const int height, const lbm::BoundaryCondition boundary)
    {
        return height + (boundary == lbm::BoundaryCondition::OnGrid ? 1 : 2);
    }
    // Obstacle-flow case test
    bool is_obstacle_case(const lbm::CaseType case_type)
    {
        return case_type == lbm::CaseType::Cylinder || case_type == lbm::CaseType::Airfoil;
    }
    // Cylinder stability diagnostic message
    std::string stability_message(const lbm::Config &cfg)
    {
        const int length = cfg.case_type == lbm::CaseType::Airfoil ? lbm::airfoil_chord(cfg) : 2 * lbm::cylinder_radius(cfg);
        std::ostringstream message;
        message << lbm::case_name(cfg.case_type) << " setup is outside the safe LBM range: Re= " << cfg.reynolds_number
                << ", grid= " << cfg.ny << ", length= " << length << ", inlet_ux= " << lbm::cylinder_inlet_ux(cfg)
                << ", nu= " << lbm::viscosity(cfg) << ", tau= " << lbm::relaxation_time(cfg) << ". ";
        return message.str();
    }
    // Parse wall boundary condition
    lbm::BoundaryCondition parse_boundary(const std::string &value)
    {
        if (value == "on-grid" || value == "ongrid")
        {
            return lbm::BoundaryCondition::OnGrid;
        }
        if (value == "mid-grid" || value == "midgrid" || value == "half-way" || value == "halfway")
        {
            return lbm::BoundaryCondition::MidGrid;
        }
        throw std::invalid_argument("invalid boundary condition: " + value + " (expected on-grid or mid-grid)");
    }
    // Parse simulation case
    lbm::CaseType parse_case(const std::string &value)
    {
        if (value == "poiseuille")
        {
            return lbm::CaseType::Poiseuille;
        }
        if (value == "cylinder")
        {
            return lbm::CaseType::Cylinder;
        }
<<<<<<< HEAD
        if (value == "convection")
        {
            return lbm::CaseType::Convection;
        }
        throw std::invalid_argument("invalid case: " + value + " (expected poiseuille, cylinder, or convection)");
=======
        if (value == "airfoil")
        {
            return lbm::CaseType::Airfoil;
        }
        throw std::invalid_argument("invalid case: " + value + " (expected poiseuille, cylinder, or airfoil)");
>>>>>>> 58006a0acf6f26ac49752afa8b4c9d4567b46eb1
    }
    // Command-line help
    void print_usage()
    {
        std::cout << "Usage: lbm_solver [--grid N|NxM] [--nx N] [--ny N] [--steps N] [--nu NU]\n"
<<<<<<< HEAD
                  << "                  [--case poiseuille|cylinder|convection]\n"
=======
                  << "                  [--case poiseuille|cylinder|airfoil]\n"
>>>>>>> 58006a0acf6f26ac49752afa8b4c9d4567b46eb1
                  << "                  [--force-x F] [--boundary on-grid|mid-grid]\n"
                  << "                  [--re RE] [--inlet-ux U] [--outlet-rho RHO]\n"
                  << "                  [--kappa K] [--gravity G] [--beta B] [--th Th] [--tc Tc]\n"
                  << "                  [--cylinder-x X] [--cylinder-y Y] [--cylinder-radius R]\n"
                  << "                  [--airfoil-x X] [--airfoil-y Y] [--airfoil-chord C] [--airfoil-angle A]\n"
                  << "                  [--report-interval N] [--output path]\n"
                  << "                  [--snapshot-interval N] [--snapshot-dir path]\n"
                  << "                  [--threads N]\n";
    }
    // Parse command-line arguments into Config
    lbm::Config parse_args(const int argc, char **argv)
    {
        lbm::Config cfg;
        bool square_grid_was_set = false;
        bool rectangular_grid_was_set = false;
        bool explicit_ny_was_set = false;
        bool explicit_nu_was_set = false;
        bool force_x_was_set = false;
        bool output_was_set = false;
        bool snapshot_dir_was_set = false;
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
                square_grid_was_set = parse_grid(require_value(arg), cfg);
                rectangular_grid_was_set = !square_grid_was_set;
            }
            else if (arg == "--case")
            {
                cfg.case_type = parse_case(require_value(arg));
            }
            else if (arg == "--nx")
            {
                cfg.nx = parse_int(require_value(arg), arg);
            }
            else if (arg == "--ny")
            {
                cfg.ny = parse_int(require_value(arg), arg);
                explicit_ny_was_set = true;
            }
            else if (arg == "--steps")
            {
                cfg.steps = parse_int(require_value(arg), arg);
            }
            else if (arg == "--nu")
            {
                cfg.nu = parse_double(require_value(arg), arg);
                explicit_nu_was_set = true;
            }
            else if (arg == "--tau")
            {
                const double tau = parse_double(require_value(arg), arg);
                if (tau <= 0.5)
                {
                    throw std::invalid_argument("tau must be greater than 0.5 for positive viscosity");
                }
                cfg.nu = lbm::cs2 * (tau - 0.5);
                explicit_nu_was_set = true;
            }
            else if (arg == "--force-x")
            {
                cfg.force_x = parse_double(require_value(arg), arg);
                force_x_was_set = true;
            }
            else if (arg == "--re")
            {
                cfg.reynolds_number = parse_double(require_value(arg), arg);
            }
            else if (arg == "--inlet-ux")
            {
                cfg.inlet_ux = parse_double(require_value(arg), arg);
            }
            else if (arg == "--outlet-rho")
            {
                cfg.outlet_rho = parse_double(require_value(arg), arg);
            }
            else if (arg == "--cylinder-x")
            {
                cfg.cylinder_x = parse_int(require_value(arg), arg);
            }
            else if (arg == "--cylinder-y")
            {
                cfg.cylinder_y = parse_int(require_value(arg), arg);
            }
            else if (arg == "--cylinder-radius")
            {
                cfg.cylinder_radius = parse_int(require_value(arg), arg);
            }
            else if (arg == "--airfoil-x")
            {
                cfg.airfoil_x = parse_int(require_value(arg), arg);
            }
            else if (arg == "--airfoil-y")
            {
                cfg.airfoil_y = parse_int(require_value(arg), arg);
            }
            else if (arg == "--airfoil-chord")
            {
                cfg.airfoil_chord = parse_int(require_value(arg), arg);
            }
            else if (arg == "--airfoil-angle")
            {
                cfg.airfoil_angle = parse_double(require_value(arg), arg);
            }
            else if (arg == "--boundary")
            {
                cfg.boundary = parse_boundary(require_value(arg));
            }
            // --- 對流相關參數 ---
            else if (arg == "--kappa")
            {
                cfg.kappa = parse_double(require_value(arg), arg);
            }
            else if (arg == "--gravity")
            {
                cfg.gravity = parse_double(require_value(arg), arg);
            }
            else if (arg == "--beta")
            {
                cfg.beta = parse_double(require_value(arg), arg);
            }
            else if (arg == "--th")
            {
                cfg.Th = parse_double(require_value(arg), arg);
            }
            else if (arg == "--tc")
            {
                cfg.Tc = parse_double(require_value(arg), arg);
            }
            // --------------------
            else if (arg == "--report-interval")
            {
                cfg.report_interval = parse_int(require_value(arg), arg);
            }
            else if (arg == "--output")
            {
                cfg.output = require_value(arg);
                output_was_set = true;
            }
            else if (arg == "--snapshot-interval")
            {
                cfg.snapshot_interval = parse_int(require_value(arg), arg);
            }
            else if (arg == "--threads")
            {
                cfg.threads = parse_int(require_value(arg), arg);
            }
            else if (arg == "--snapshot-dir")
            {
                cfg.snapshot_dir = require_value(arg);
                snapshot_dir_was_set = true;
            }
            else if (arg == "--help")
            {
                print_usage();
                std::exit(0);
            }
            else
            {
                throw std::invalid_argument("unknown option: " + arg);
            }
        }

<<<<<<< HEAD
        // 設定各個 case 的預設輸出目錄與檔案
        if (cfg.case_type == lbm::CaseType::Cylinder)
=======
        if (is_obstacle_case(cfg.case_type))
>>>>>>> 58006a0acf6f26ac49752afa8b4c9d4567b46eb1
        {
            if (!output_was_set)
            {
                cfg.output = cfg.case_type == lbm::CaseType::Airfoil ? "test/airfoil_field.csv" : "test/cylinder_field.csv";
            }
            if (!snapshot_dir_was_set)
            {
                cfg.snapshot_dir = cfg.case_type == lbm::CaseType::Airfoil ? "test/airfoil_snapshots" : "test/cylinder_snapshots";
            }
        }
        else if (cfg.case_type == lbm::CaseType::Convection)
        {
            if (!output_was_set)
            {
                cfg.output = "test/convection_field.csv";
            }
            if (!snapshot_dir_was_set)
            {
                cfg.snapshot_dir = "test/convection_snapshots";
            }
        }
        else
        {
            // Poiseuille 預設路徑 (假設原本如果沒設定會交給內部處理，或統一設在這裡)
            if (!snapshot_dir_was_set)
            {
                cfg.snapshot_dir = "test/poiseuille_snapshots";
            }
        }

<<<<<<< HEAD
        // 網格調整邏輯
        if (square_grid_was_set && cfg.case_type == lbm::CaseType::Cylinder && !explicit_ny_was_set)
=======
        if (square_grid_was_set && is_obstacle_case(cfg.case_type) && !explicit_ny_was_set)
>>>>>>> 58006a0acf6f26ac49752afa8b4c9d4567b46eb1
        {
            const int grid_size = cfg.nx;
            cfg.nx = 4 * grid_size;
            cfg.ny = grid_size;
        }
        if (cfg.case_type == lbm::CaseType::Poiseuille && !rectangular_grid_was_set && !explicit_ny_was_set)
        {
            const int channel_height = square_grid_was_set ? cfg.nx : lbm::default_grid_size;
            cfg.ny = poiseuille_y_nodes(channel_height, cfg.boundary);
        }
        if (cfg.case_type == lbm::CaseType::Poiseuille && !force_x_was_set)
        {
            cfg.force_x = lbm::default_pressure_drop / static_cast<double>(cfg.nx);
        }
        if (is_obstacle_case(cfg.case_type))
        {
            if (explicit_nu_was_set)
            {
                throw std::invalid_argument("Obstacle-flow viscosity is inferred from --re and --inlet-ux; remove --nu or --tau");
            }
            if (cfg.reynolds_number <= 0.0)
            {
                throw std::invalid_argument("Obstacle-flow Reynolds number must be positive");
            }
            if (cfg.inlet_ux <= 0.0)
            {
                throw std::invalid_argument("Obstacle-flow inlet velocity must be positive");
            }
            cfg.nu = cfg.case_type == lbm::CaseType::Airfoil ? lbm::airfoil_viscosity_from_re(cfg) : lbm::cylinder_viscosity_from_re(cfg);
        }

        // 基本防呆檢查
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
            throw std::invalid_argument("Snapshot interval must be non-negative");
        }
        if (cfg.threads < 0)
        {
            throw std::invalid_argument("Thread count must be non-negative");
        }
        if (cfg.nu <= 0.0)
        {
            throw std::invalid_argument("Kinematic viscosity nu must be positive");
        }
<<<<<<< HEAD

        // Cylinder 防呆檢查
        if (cfg.case_type == lbm::CaseType::Cylinder)
=======
        if (is_obstacle_case(cfg.case_type))
>>>>>>> 58006a0acf6f26ac49752afa8b4c9d4567b46eb1
        {
            if (cfg.nx < 10 || cfg.ny < 10)
            {
                throw std::invalid_argument("Obstacle-flow case needs a larger domain");
            }
            if (cfg.inlet_ux < lbm::min_cylinder_inlet_ux || cfg.inlet_ux > lbm::max_cylinder_inlet_ux)
            {
                std::ostringstream message;
                message << stability_message(cfg) << "Use " << lbm::min_cylinder_inlet_ux << " <= --inlet-ux <= "
                        << lbm::max_cylinder_inlet_ux << ".";
                throw std::invalid_argument(message.str());
            }
            if (lbm::relaxation_time(cfg) < lbm::min_cylinder_tau)
            {
                std::ostringstream message;
                message << stability_message(cfg) << "tau is too close to 0.5. Use a larger grid, lower Re, or larger "
                        << "--inlet-ux within the safe velocity range.";
                throw std::invalid_argument(message.str());
            }
            if (lbm::relaxation_time(cfg) > lbm::max_cylinder_tau)
            {
                std::ostringstream message;
                message << stability_message(cfg) << "tau is too large. Use a smaller grid, larger Re, or smaller "
                        << "--inlet-ux within the safe velocity range.";
                throw std::invalid_argument(message.str());
            }
            if (cfg.outlet_rho <= 0.0)
            {
                throw std::invalid_argument("Outlet density must be positive");
            }
        }
        
        // Convection 防呆檢查
        if (cfg.case_type == lbm::CaseType::Convection)
        {
            if (cfg.kappa <= 0.0)
            {
                throw std::invalid_argument("Thermal diffusivity kappa must be positive");
            }
        }

        return cfg;
    }
} // namespace

int main(int argc, char **argv)
{
    try
    {
        const lbm::Config cfg = parse_args(argc, argv);
        lbm::configure_threads(cfg);

        switch (cfg.case_type)
        {
        case lbm::CaseType::Poiseuille:
            lbm::run_poiseuille(cfg);
            break;
        case lbm::CaseType::Cylinder:
            lbm::run_cylinder(cfg);
            break;
<<<<<<< HEAD
        case lbm::CaseType::Convection:
            lbm::run_convection(cfg);
=======
        case lbm::CaseType::Airfoil:
            lbm::run_airfoil(cfg);
>>>>>>> 58006a0acf6f26ac49752afa8b4c9d4567b46eb1
            break;
        }
    }
    catch (const std::exception &error)
    {
        std::cerr << "error: " << error.what() << '\n';
        return 1;
    }
    return 0;
}