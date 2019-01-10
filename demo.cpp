/*
int main()
{
	double a = 3.0, b = 5.0, c = 1.0, d = 2;        // 真实参数值
	int N = 100;                          // 数据点
	double w_sigma = 1.0 / 20;                 // 噪声Sigma值
	cv::RNG rng;                        // OpenCV随机数产生器
	double abcd[4] = { 1,2, 1,1 };            // abc参数的估计值

	vector<double> x_data, y_data;      // 数据

	cout << "generating data: " << endl;
	for (int i = 0; i < N; i++)
	{
		double x = i + 0.2;
		x_data.push_back(x);
		y_data.push_back(a * pow(1 + exp(b - c * x), -1.0 / d) + rng.gaussian(w_sigma));
		cout << x_data[i] << " " << y_data[i] << endl;
	}

	typedef g2o::BlockSolver<g2o::BlockSolverTraits<4, 1>> Block;
	Block::LinearSolverType* linearSolver = new g2o::LinearSolverDense<Block::PoseMatrixType>();
	Block* solver_ptr = new Block(linearSolver);

	g2o::OptimizationAlgorithmLevenberg* solver = new g2o::OptimizationAlgorithmLevenberg(solver_ptr);
	g2o::SparseOptimizer optimizer;
	optimizer.setAlgorithm(solver);
	optimizer.setVerbose(true);

	VertexCurveFitting* v = new VertexCurveFitting();
	v->setEstimate(Eigen::Matrix<double, 4, 1>{1,1,2,6});
	v->setId(0);
	optimizer.addVertex(v);

	for (int i = 0; i < N; i++)
	{
		EdgeCurveFitting* edge = new EdgeCurveFitting(x_data[i]);
		edge->setId(i+1);

		edge->setVertex(0, v);
		edge->setMeasurement(y_data[i]);
		edge->setInformation(Eigen::Matrix<double, 1, 1>::Identity() * 1 / (w_sigma*w_sigma));
		optimizer.addEdge(edge);
	}

	// 执行优化
	cout << "start optimization" << endl;
	chrono::steady_clock::time_point t1 = chrono::steady_clock::now();
	optimizer.initializeOptimization();
	optimizer.optimize(100);
	chrono::steady_clock::time_point t2 = chrono::steady_clock::now();
	chrono::duration<double> time_used = chrono::duration_cast<chrono::duration<double>>(t2 - t1);
	cout << "solve time cost = " << time_used.count() << " seconds. " << endl;

	// 输出优化值
	Eigen::Vector4d abc_estimate = v->estimate();
	cout << "estimated model: " << abc_estimate.transpose() << endl;

	return 0;
}

int main()
{
	double a = 3.0, b = 5.0, c = 1.0, d = 2;        // 真实参数值
	int N = 100;                          // 数据点
	double w_sigma = 1.0/20;                 // 噪声Sigma值
	cv::RNG rng;                        // OpenCV随机数产生器
	double abcd[4] = { 1,2, 1,1 };            // abc参数的估计值

	vector<double> x_data, y_data;      // 数据

	cout << "generating data: " << endl;
	for (int i = 0; i < N; i++)
	{
		double x = i+0.2 ;
		x_data.push_back(x);
		y_data.push_back(a * pow(1 + exp(b - c * x), -1.0 / d) + rng.gaussian(w_sigma));
		cout << x_data[i] << " " << y_data[i] << endl;
	}
	ceres::Problem problem;

	for (int i = 0; i < N; i++)
	{
		Rat43Analytic* factor = new Rat43Analytic(x_data[i], y_data[i]);
		problem.AddResidualBlock(factor, NULL, abcd);
	}
	ceres::Solver::Options options;
	options.max_num_iterations = 1000;
	options.linear_solver_type = ceres::DENSE_QR;
	ceres::Solver::Summary summary;
	Solve(options, &problem, &summary);

	cout << "a=:"<<abcd[0]<< endl << "b=:"<< abcd[1]<< endl << "c=:" << abcd[2]<< 
		endl << "d=:" << abcd[3] << endl;
	return 0;
}
*/
