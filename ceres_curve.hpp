#pragma once
#include <Eigen/Core>
#include <ceres/ceres.h>
#include <ceres/cost_function.h>

class Rat43Analytic :public ceres::SizedCostFunction<1, 4>
{
public:
	Rat43Analytic(const double x, const double y) :x_(x), y_(y) {}
	virtual ~Rat43Analytic() {}

	virtual bool Evaluate(double const* const* parameters, double* residuals,
		double** jacobians) const
	{
		//const double a = parameters[0][0];
		//const double b = parameters[0][1];
		//const double c = parameters[0][2];
		//const double d = parameters[0][3];
		const double b1 = parameters[0][0];
		const double b2 = parameters[0][1];
		const double b3 = parameters[0][2];
		const double b4 = parameters[0][3];


		//residuals[0] =  std::exp(a* _x + b * sin(_x) + c) * d - _y ;
		//residuals[0] = b1 * pow(1 + exp(b2 - b3 * x_), -1.0 / b4) - y_;

		const double t1 = exp(b2 - b3 * x_);
		const double t2 = 1 + t1;
		const double t3 = pow(t2, -1.0 / b4);
		residuals[0] = b1 * t3 - y_;

		if (!jacobians) return true;
		double* jacobian = jacobians[0];
		if (!jacobian) return true;
		//jacobian[0] = _x * std::exp(a* _x + b * sin(_x) + c) * d;
		//jacobian[1] = sin(_x) * std::exp(a* _x + b * sin(_x) + c) * d;
		//jacobian[2] = std::exp(a* _x + b * sin(_x) + c) * d;
		//jacobian[3] = std::exp(a* _x + b * sin(_x) + c);
		//// Compute the Jacobian if asked for.
		//return true;

		const double t4 = pow(t2, -1.0 / b4 - 1);
		jacobian[0] = t3;
		jacobian[1] = -b1 * t1 * t4 / b4;
		jacobian[2] = -x_ * jacobian[1];
		jacobian[3] = b1 * log(t2) * t3 / (b4 * b4);
		return true;
	}
private:
	const double x_, y_;
};
