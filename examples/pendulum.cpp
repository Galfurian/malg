/// @file pendulum.cpp
/// @author Enrico Fraccaroli (enry.frak@gmail.com)

#include <malg/control/control.hpp>
#include <malg/io.hpp>

#include <stopwatch/stopwatch.hpp>

#include <solver/stepper/stepper_adaptive.hpp>
#include <solver/stepper/stepper_euler.hpp>
#include <solver/stepper/stepper_rk4.hpp>
#include <solver/detail/observer.hpp>
#include <solver/solver.hpp>

#ifdef MALG_ENABLE_PLOT
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

    Parameter(Variable _m = 3.0,
              Variable _l = 1.19,
              Variable _b = 0.75,
              Variable _g = 9.81)
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

    Model(Parameter parameter = Parameter())
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
    inline void operator()(const State &x, State &dxdt, Time t) noexcept
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

} // namespace pendulum

template <std::size_t DECIMATION = 0>
struct ObserverSave : public solver::detail::DecimationObserver<DECIMATION> {
    std::vector<pendulum::Variable> time, angle, velocity;
    ObserverSave() = default;
    inline void operator()(const pendulum::State &x, const pendulum::Time &t) noexcept
    {
        if (this->observe()) {
            time.emplace_back(t);
            angle.emplace_back(x[0]);
            velocity.emplace_back(x[1]);
        }
    }
};

int main(int, char *[])
{
    using namespace pendulum;
    // Instantiate the model.
    Model system;
    // Runtime state.
    State x, dx;
    // Initial and runtime states.
    const State x0{ .0, .0 };
    // Simulation parameters.
    const Time time_start = 0.0;
    const Time time_end   = 40.0;
    const Time time_delta = 0.01;
    // Setup the adaptive solver.
    using FStepper        = solver::stepper_rk4<State, Time>;
    const auto Iterations = 2;
    const auto Error      = solver::ErrorFormula::Mixed;
    using AStepper        = solver::stepper_adaptive<FStepper, Iterations, Error>;
    // Instantiate the solvers.
    FStepper fstepper;
    AStepper astepper;
    astepper.set_tollerance(1e-03);
    astepper.set_min_delta(0.0001);
    astepper.set_max_delta(0.1);
    // Instantiate the observers.
#ifdef MALG_ENABLE_PLOT
    ObserverSave<0> fobserver;
    ObserverSave<0> aobserver;
#elif 1
    solver::detail::ObserverPrint<0> fobserver;
    solver::detail::ObserverPrint<0> aobserver;
#endif

    // Instantiate the stopwatch.
    stopwatch::Stopwatch sw;
    std::cout << std::fixed;
    std::cout << "Simulating...\n";
    
    // Start the simulation.
    sw.start();
    solver::integrate_fixed(fstepper, fobserver, system, x = x0, time_start, time_end, time_delta);
    sw.round();
    solver::integrate_adaptive(astepper, aobserver, system, x = x0, time_start, time_end, time_delta);
    sw.round();

    std::cout << "\n";
    std::cout << "Integration steps and elapsed times:\n";
    std::cout << "    Fixed stepper computed    " << std::setw(12) << fstepper.steps() << " steps, for a total of " << sw[0] << "\n";
    std::cout << "    Adaptive stepper computed " << std::setw(12) << astepper.steps() << " steps, for a total of " << sw[1] << "\n";

#ifdef MALG_ENABLE_PLOT
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
