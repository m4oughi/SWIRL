#include "../src/WindField.h"
#include <iostream>

int main() {
    // Set up test parameters for the wind fields
    Parameters parameters;
    parameters["Um"] = 10.0;
    parameters["rm"] = 5.0;
    parameters["zm"] = 2.0;
    parameters["S"] = 0.5;
    parameters["gamma"] = 1.5;
    parameters["rho0"] = 1.225;
    parameters["initial_center"] = std::vector<double>{0.0, 0.0, 0.0};
    parameters["initial_velocity"] = std::vector<double>{1.0, 1.0, 0.0};

    // Test BakerSterlingVortex model
    WindField* windField = new_WindField("BakerSterlingVortex", parameters);
    if (windField) {
        const int num_points = 3;
        double x[num_points] = {1.0, 2.0, 3.0};
        double y[num_points] = {1.0, 2.0, 3.0};
        double z[num_points] = {0.5, 1.0, 1.5};
        double vx[num_points], vy[num_points], vz[num_points], rhof[num_points];

        windField->get_fluid_velocity_and_density(num_points, 1.0, x, y, z, vx, vy, vz, rhof);

        std::cout << "BakerSterlingVortex velocities and densities at time 1.0:\n";
        for (int i = 0; i < num_points; i++) {
            std::cout << "Point " << i + 1 << ": vx = " << vx[i] << ", vy = " << vy[i]
                      << ", vz = " << vz[i] << ", density = " << rhof[i] << '\n';
        }
        delete windField;
    } else {
        std::cerr << "Failed to create BakerSterlingVortex model.\n";
    }

    // Add additional parameters specific to RankineVortex
    parameters["rc"] = 3.0;
    parameters["E"] = 0.8;

    // Test RankineVortex model
    windField = new_WindField("RankineVortex", parameters);
    if (windField) {
        const int num_points = 3;
        double x[num_points] = {1.0, 2.0, 3.0};
        double y[num_points] = {1.0, 2.0, 3.0};
        double z[num_points] = {0.5, 1.0, 1.5};
        double vx[num_points], vy[num_points], vz[num_points], rhof[num_points];

        windField->get_fluid_velocity_and_density(num_points, 1.0, x, y, z, vx, vy, vz, rhof);

        std::cout << "RankineVortex velocities and densities at time 1.0:\n";
        for (int i = 0; i < num_points; i++) {
            std::cout << "Point " << i + 1 << ": vx = " << vx[i] << ", vy = " << vy[i]
                      << ", vz = " << vz[i] << ", density = " << rhof[i] << '\n';
        }
        delete windField;
    } else {
        std::cerr << "Failed to create RankineVortex model.\n";
    }

    return 0;
}
