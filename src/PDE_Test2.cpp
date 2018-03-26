#include <PDE.hpp>
#include <array>
#include <fstream>

int main()
{
    // Using directives
    using real = double;
    using solver_internal = real;
    using vec2 = PDE::StateVector<real, real>;
    using body = PDE::StateVector<vec2, vec2>;
    using state = PDE::StateVector<body, body>;
    using solver = PDE::RK4::Solver<solver_internal, state>;

    // Simulation params
    real M = 100.0,
         m = 1.0,
         G = 1.0;
    solver_internal dt = 0.05,
                    max_t = 200.0;

    solver rk4;

    //                    pos    ,     vel
    rk4.lhs() = state{{{0.0, 0.0}, {0.0, 0.0}},
                      {{1.0, 1.0}, {0.0, 0.05}}};

    rk4.equation() = [=](state& result, const state& rhs)
    {
        const auto& central_body = rhs.get<0>();
        const auto& outer_body = rhs.get<1>();

        auto dist = [](const body& body1, const body& body2)
        {
            auto body1_pos = body1.get<0>();
            auto body2_pos = body2.get<0>();

            return std::sqrt(std::pow(body1_pos.get<0>() - body2_pos.get<0>(), 2) +
                             std::pow(body1_pos.get<1>() - body2_pos.get<1>(), 2));
        };

        result = PDE::make_equation(body{m * G * (outer_body.get<0>() - central_body.get<0>()) / std::pow(dist(central_body, outer_body), 3), },
                                    body{M * G * (central_body.get<0>() - outer_body.get<0>()) / std::pow(dist(central_body, outer_body), 3), });
    };

    std::ofstream out("newton.dat");

    for(solver_internal t = 0.0 ; t < max_t ; t += dt)
    {
        out << t << '\t' << rk4.lhs().get<0>().get<0>().get<0>() << '\t' << rk4.lhs().get<0>().get<0>().get<1>() << 
                    '\t' << rk4.lhs().get<1>().get<0>().get<0>() << '\t' << rk4.lhs().get<1>().get<0>().get<1>() << std::endl;

        rk4.iterate(dt);
    }

    return 0;
}