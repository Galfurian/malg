/// @file pendulum.cpp
/// @author Enrico Fraccaroli (enry.frak@gmail.com)

#include <malg/control/control.hpp>
#include <malg/io.hpp>

#include <stopwatch/stopwatch.hpp>

#include <chainsaw/detail/observer.hpp>
#include <chainsaw/solver.hpp>
#include <chainsaw/stepper/stepper_adaptive.hpp>
#include <chainsaw/stepper/stepper_euler.hpp>
#include <chainsaw/stepper/stepper_rk4.hpp>

#ifdef ENABLE_PLOT
#include <matplot/matplot.h>
#endif

namespace pendulum
{

/// The time is a continuous time value.
using Time = double;
/// State variables are declared as floating point.
using Variable = double;
/// Matrix type
using Matrix = malg::Matrix<Variable>;
/// Vector type.
using Vector = malg::Vector<Variable>;
/// @brief State of the system.
///     x[0] : Angle.
///     x[1] : Velocity.
using State = Vector;

/// @brief Parameters of our model.
struct Parameter {
    /// @brief Mass of the rod [kg].
    const Variable m;
    /// @brief Lenght of the rod [m].
    const Variable l;
    /// @brief Rotational damping [N.m].
    const Variable b;
    /// @brief Gravitational force [N].
    const Variable g;

    explicit Parameter()
        : Parameter(3.0, 1.19, 0.75, 9.81)
    {
        // Nothing to do.
    }

    explicit Parameter(Variable _m, Variable _l, Variable _b, Variable _g)
        : m(_m),
          l(_l),
          b(_b),
          g(_g)
    {
        // Nothing to do.
    }
};

struct Model : public Parameter {
    /// Continous time state-space model.
    malg::control::StateSpace<Variable> sys;
    /// Input matrix.
    Vector u;

    explicit Model(Parameter parameter)
        : Parameter(parameter),
          sys(),
          u(1, .0)
    {
        sys.A = {
            { 0., 1. },
            { -g / l, -b / (m * l * l) }
        };
        sys.B = {
            { 0. },
            { 1. / (m * l * l) }
        };
        sys.C = {
            { 0., 1. }
        };
        sys.D = {
            { 0. }
        };
    }

    /// @brief DC motor behaviour.
    /// @param x the current state.
    /// @param dxdt the final state.
    /// @param t the current time.
    inline void operator()(const State &x, State &dxdt, Time t)
    {
#if 1
        u[0] = (t < 5) ? 1 : 0;
#else
        u[0] = .0;
#endif
        // Equation system.
        dxdt = malg::dot(sys.A, x) + malg::dot(sys.B, u);
    }
};

template <std::size_t DECIMATION = 0>
struct ObserverSave : public chainsaw::detail::ObserverDecimate<State, Time, DECIMATION> {
    inline void operator()(const State &x, const Time &t) noexcept override
    {
        if (this->observe()) {
            time.emplace_back(t);
            angle.emplace_back(x[0]);
            velocity.emplace_back(x[1]);
        }
    }

    std::vector<Variable> time, angle, velocity;
};

} // namespace pendulum

int main(int, char *[])
{
    using namespace pendulum;
    // Instantiate the parameters.
    Parameter parameters;
    // Instantiate the model.
    Model system(parameters);
    // Runtime state.
    State x, dx;
    // Initial and runtime states.
    const State x0{ .0, .0 };
    // Simulation parameters.
    const Time time_start = 0.0;
    const Time time_end   = 40.0;
    const Time time_delta = 0.01;
    // Setup the adaptive solver.
    using FStepper        = chainsaw::stepper_rk4<State, Time>;
    const auto Iterations = 2;
    const auto Error      = chainsaw::ErrorFormula::Mixed;
    using AStepper        = chainsaw::stepper_adaptive<FStepper, Iterations, Error>;
    // Instantiate the solvers.
    FStepper fstepper;
    AStepper astepper;
    astepper.set_tollerance(1e-03);
    astepper.set_min_delta(0.0001);
    astepper.set_max_delta(0.1);

    // Instantiate the observers.
#ifdef ENABLE_PLOT
    using Observer = ObserverSave<0>;
#else
    using Observer = chainsaw::detail::ObserverPrint<State, Time, 0>;
#endif

    Observer fobserver;
    Observer aobserver;

    // Instantiate the stopwatch.
    stopwatch::Stopwatch sw;
    std::cout << std::fixed;
    std::cout << "Simulating...\n";

    // Start the simulation.
    sw.start();
    chainsaw::integrate_fixed(fstepper, fobserver, system, x = x0, time_start, time_end, time_delta);
    sw.round();
    chainsaw::integrate_adaptive(astepper, aobserver, system, x = x0, time_start, time_end, time_delta);
    sw.round();

    std::cout << "\n";
    std::cout << "Integration steps and elapsed times:\n";
    std::cout << "    Fixed stepper computed    " << std::setw(12) << fstepper.steps() << " steps, for a total of " << sw[0] << "\n";
    std::cout << "    Adaptive stepper computed " << std::setw(12) << astepper.steps() << " steps, for a total of " << sw[1] << "\n";

#ifdef ENABLE_PLOT
    auto figure = matplot::figure(true);
    matplot::grid(matplot::on);
    matplot::hold(matplot::on);
    matplot::plot(fobserver.time, fobserver.angle)->line_width(2).display_name("[F] Angle A (rad)");
    matplot::plot(aobserver.time, aobserver.angle)->line_width(2).display_name("[A] Angle A (rad)");
    matplot::plot(fobserver.time, fobserver.velocity, "--")->line_width(1).display_name("[F] Angular Speed A (rad/s)");
    matplot::plot(aobserver.time, aobserver.velocity, "--")->line_width(1).display_name("[A] Angular Speed A (rad/s)");
    matplot::xlabel("Time (s)");
    matplot::legend(matplot::on)->location(matplot::legend::general_alignment::top);
    matplot::show();
#endif
    return 0;
}
