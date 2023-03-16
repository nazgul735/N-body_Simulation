#include <iostream>
#include <vector>
#include <cmath>

const double G = 6.6743e-11; // gravitational constant

class Body {
public:
    double mass;
    double x, y, z;
    double vx, vy, vz;
};

class NBodySimulator {
public:
    std::vector<Body> bodies;
    void simulate(double dt, double t_max);
private:
    void update(double dt);
    void computeForces();
    double distance(Body b1, Body b2);
    double force(double m1, double m2, double r);
    double acceleration(double f, double m);
};

void NBodySimulator::simulate(double dt, double t_max) {
    double t = 0;
    while (t < t_max) {
        computeForces();
        update(dt);
        t += dt;
    }
}

void NBodySimulator::computeForces() {
    for (int i = 0; i < bodies.size(); i++) {
        Body& b1 = bodies[i];
        for (int j = i+1; j < bodies.size(); j++) {
            Body& b2 = bodies[j];
            double r = distance(b1, b2);
            double f = force(b1.mass, b2.mass, r);
            double a1 = acceleration(f, b1.mass);
            double a2 = acceleration(f, b2.mass);
            double dx = b2.x - b1.x;
            double dy = b2.y - b1.y;
            double dz = b2.z - b1.z;
            b1.vx += a1 * dx / r;
            b1.vy += a1 * dy / r;
            b1.vz += a1 * dz / r;
            b2.vx -= a2 * dx / r;
            b2.vy -= a2 * dy / r;
            b2.vz -= a2 * dz / r;
        }
    }
}

void NBodySimulator::update(double dt) {
    for (auto& body : bodies) {
        body.x += body.vx * dt;
        body.y += body.vy * dt;
        body.z += body.vz * dt;
    }
}

double NBodySimulator::distance(Body b1, Body b2) {
    double dx = b2.x - b1.x;
    double dy = b2.y - b1.y;
    double dz = b2.z - b1.z;
    return std::sqrt(dx*dx + dy*dy + dz*dz);
}

double NBodySimulator::force(double m1, double m2, double r) {
    return G * m1 * m2 / (r*r);
}

double NBodySimulator::acceleration(double f, double m) {
    return f / m;
}

int main() {
    // Create a simulator and some bodies
    NBodySimulator simulator;
    Body sun = {1.989e30, 0, 0, 0, 0, 0};
    Body earth = {5.972e24, 149.6e9, 0, 0, 29.78e3, 0};
    Body moon = {7.342e22, 149.6e9 + 384.4e6, 0, 0, 29.78e3 + 1.02e3, 0};
    simulator.bodies.push_back(sun);
    simulator.bodies.push_back(earth);
    simulator.bodies.push_back(moon);
    
    //