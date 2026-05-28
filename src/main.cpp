#include "lbm.hpp"

#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>

namespace
{
    //
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
    //
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
    //
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
    //
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
    //
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
        throw std::invalid_argument("invalid case: " + value + " (expected poiseuille or cylinder)");
    }
    //
    void print_usage()
    {
        std::cout << "Usage: lbm_solver [--grid N|NxM] [--nx N] [--ny N] [--steps N] [--nu NU]\n"
                  << "                  [--case poiseuille|cylinder]\n"
                  << "                  [--force-x F] [--boundary on-grid|mid-grid]\n"
                  << "                  [--re RE] [--outlet-rho RHO]\n"
                  << "                  [--cylinder-x X] [--cylinder-y Y] [--cylinder-radius R]\n"
                  << "                  [--report-interval N] [--output path]\n"
                  << "                  [--snapshot-interval N] [--snapshot-dir path]\n"
                  << "                  [--threads N]\n";
    }
    //
    lbm::Config parse_args(const int argc, char **argv)
    {
        lbm::Config cfg;
        bool square_grid_was_set = false;
        bool explicit_ny_was_set = false;
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
            }
            else if (arg == "--tau")
            {
                const double tau = parse_double(require_value(arg), arg);
                if (tau <= 0.5)
                {
                    throw std::invalid_argument("tau must be greater than 0.5 for positive viscosity");
                }
                cfg.nu = lbm::cs2 * (tau - 0.5);
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

        if (cfg.case_type == lbm::CaseType::Cylinder)
        {
            if (!output_was_set)
            {
                cfg.output = "test/cylinder_field.csv";
            }
            if (!snapshot_dir_was_set)
            {
                cfg.snapshot_dir = "test/cylinder_snapshots";
            }
        }
        else if (!snapshot_dir_was_set)
        {
            cfg.snapshot_dir = "test/poiseuille_snapshots";
        }

        if (square_grid_was_set && cfg.case_type == lbm::CaseType::Cylinder && !explicit_ny_was_set)
        {
            const int grid_size = cfg.nx;
            cfg.nx = 4 * grid_size;
            cfg.ny = grid_size;
        }
        if (cfg.case_type == lbm::CaseType::Poiseuille && !force_x_was_set)
        {
            cfg.force_x = lbm::default_pressure_drop / static_cast<double>(cfg.nx);
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
        if (cfg.case_type == lbm::CaseType::Cylinder)
        {
            if (cfg.nx < 10 || cfg.ny < 10)
            {
                throw std::invalid_argument("Cylinder case needs a larger domain");
            }
            if (cfg.reynolds_number <= 0.0)
            {
                throw std::invalid_argument("Cylinder Reynolds number must be positive");
            }
            if (cfg.outlet_rho <= 0.0)
            {
                throw std::invalid_argument("Outlet density must be positive");
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
        }
    }
    catch (const std::exception &error)
    {
        std::cerr << "error: " << error.what() << '\n';
        return 1;
    }
    return 0;
}
