#pragma once
#include <iostream>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <g2o/core/base_vertex.h>
#include <g2o/core/base_unary_edge.h>


class VertexCurveFitting :public g2o::BaseVertex<4, Eigen::Matrix<double, 4, 1>>
{
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	VertexCurveFitting() { setToOriginImpl(); }

	virtual void setToOriginImpl()
	{
		_estimate = Eigen::Matrix<double, 4, 1>::Zero();
	}

	virtual void oplusImpl(const double* update)
	{
		_estimate += Eigen::Matrix<double, 4, 1>(update);
	}

	virtual bool read(std::istream &in) { return true; }

	virtual bool write(std::ostream &out) const { return true; }
};


class EdgeCurveFitting : public g2o::BaseUnaryEdge<1, double, VertexCurveFitting>
{
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	EdgeCurveFitting(double x) :BaseUnaryEdge(), _x(x) {}

	void computeError()
	{
		const VertexCurveFitting* v = static_cast<const VertexCurveFitting*>(_vertices[0]);
		const Eigen::Matrix<double, 4, 1> param = v->estimate();
		double a = param[0];
		double b = param[1];
		double c = param[2];
		double d = param[3];

		_error(0, 0) = a * pow(1 + exp(b - c * _x), -1.0 / d) - _measurement;
	}

	virtual void linearizeOplus()
	{
		const VertexCurveFitting* v = static_cast<const VertexCurveFitting*>(_vertices[0]);
		const Eigen::Matrix<double, 4, 1> param = v->estimate();
		double a = param[0];
		double b = param[1];
		double c = param[2];
		double d = param[3];

		const double t1 = exp(b - c * _x);
		const double t2 = 1 + t1;
		const double t3 = pow(t2, -1.0 / d);
		const double t4 = pow(t2, -1.0 / d - 1);
		Eigen::Matrix<double, 1, 4> jacobian;
		jacobian[0] = t3;
		jacobian[1] = -a * t1 * t4 / d;
		jacobian[2] = -_x * jacobian[1];
		jacobian[3] = a * log(t2) * t3 / (d * d);

		_jacobianOplusXi = jacobian;
	}
	virtual bool read(std::istream& in) { return true; }
	virtual bool write(std::ostream& out) const { return true; }

public:
	double _x;
};