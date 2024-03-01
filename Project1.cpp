
#define _USE_MATH_DEFINES

#include <iostream>
#include <SFML/Graphics.hpp>
#include <cmath>
#include <functional>
#include <vector>

class Pendulum {
public:
    sf::Vector2f origin;
    float length;
    double theta; // angle
    double omega; // angular velocity
    double alpha; // angular acceleration
    double mass; // in kg
    double gravity = 9.81; // m/s^2

    
    Pendulum(sf::Vector2f origin, float length)
        : origin(origin), length(length), theta(M_PI/2), omega(0), alpha(0), mass(1) {};

   
    double acceleration(double angle) const {
        return (-1*gravity / length) * std::sin(angle);
    }


};

void RK4(double& y, double& dy, double dt, const std::function<double(double)>& f) {
    double k1 = dt * f(y);
    double k2 = dt * f(y + 0.5 * k1);
    double k3 = dt * f(y + 0.5 * k2);
    double k4 = dt * f(y + k3);

    y += (k1 + 2 * k2 + 2 * k3 + k4) / 6;
    dy += dt * f(y);
}

void updatePendulum(Pendulum& pendulum, double& theta, double& omega,  double dt) {
    auto omega_dot = -pendulum.gravity / pendulum.length * std::sin(theta);
    auto theta_dot = omega;
    //updates angular velocity based on its derivative 
    RK4(omega, omega_dot, dt, [&pendulum, &theta](double) -> double {
        return -pendulum.gravity / pendulum.length * std::sin(theta); 
        });

    RK4(theta, theta_dot, dt, [&pendulum, &omega](double) -> double {
        return omega; 
        });
}



int main()
{
    sf::RenderWindow window(sf::VideoMode(800, 600), "Pendulum Simulation");
    Pendulum pendulum(sf::Vector2f(400, 300), 200); 

    while (window.isOpen()) {
        sf::Event event;
        while (window.pollEvent(event)) {
            if (event.type == sf::Event::Closed)
                window.close();
        }

        // Update Pendulum
        double dt = 0.001; // Time step, 
        updatePendulum(pendulum, pendulum.theta, pendulum.omega, dt);

        // Draw
        window.clear();
        sf::CircleShape pendulumBob(10); 
        pendulumBob.setFillColor(sf::Color::Red);
        double bobX = pendulum.origin.x + pendulum.length * std::sin(pendulum.theta);
        double bobY = pendulum.origin.y + pendulum.length * std::cos(pendulum.theta);
        pendulumBob.setPosition(static_cast<float>(bobX), static_cast<float>(bobY));
        



        // Line representing the rod
        sf::VertexArray lines(sf::Lines, 2);
        lines[0].position = pendulum.origin;
        lines[1].position = sf::Vector2f(static_cast<float>(bobX), static_cast<float>(bobY));
        lines[0].color = sf::Color::White;
        lines[1].color = sf::Color::White;

        window.draw(lines);
        window.draw(pendulumBob);
        window.display();
    }



}
