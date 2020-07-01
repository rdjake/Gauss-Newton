#include "reverse.h"

inline int GaussNewton::sign(double val) {
	return (val < 0) ? -1 : 1;
}

void GaussNewton::createPath() {

	setlocale(LC_ALL, "Russian");

	if (newWay) {

		int i;
		vector<filesystem::path> dirs;
		filesystem::path parent = (filesystem::path)"C:\\Filtration";
		while (true) {
			cout << "\n";
			dirs.clear();
			for (const auto& entry : filesystem::directory_iterator(parent))
				if (entry.is_directory())
				{
					if (entry.path().filename().string() == "mesh") goto done;
					dirs.push_back(entry.path());
				}
			for (i = 0; i < dirs.size(); i++) {
				cout << to_string(i) + " - " + dirs[i].filename().string() << endl;
			}
			cout << "\n>> ";
			cin >> i;
			parent = dirs[i];
		}
	done:
		path = parent.string();
		pathHome = parent.parent_path().string();
	}

	else {
		cout << "Enter path with model from Test_Models: ";
		string p;
		cin >> p;

		pathHome = "C:\\Filtration\\Test_models";
		path = "C:\\Filtration\\Test_models\\" + p;

	}
	setlocale(LC_ALL, "C");
	pathMesh = path + "\\mesh";
	pathReverse = path + "\\reverse";
	pathOutput = path + "\\output";
	pathOutput3D = path + "\\output3D";
	pathProperties = path + "\\properties";
}

void GaussNewton::ReadData() {
	string temp;
	int tempI;
	double tempD;
	ifstream inWells(pathProperties + "\\wellsCONS.txt");
	inWells >> numOfWells;

	WellsType.resize(numOfWells);

	int tempN;
	double tempSum;
	for (int i = 0; i < numOfWells; i++) {
		inWells >> tempI;
		inWells >> tempI;
		inWells >> temp;
		inWells >> tempN;
		tempSum = 0.0;
		for (int j = 0; j < tempN; j++) {
			inWells >> tempI;

			inWells >> tempD;
			tempSum += tempD;

			inWells >> tempI;
			inWells >> tempI;
			inWells >> tempI;

		}
		if (tempSum > 0.0) {
			WellsType[i] = false;
			numOfIn++;
		}
		else {
			WellsType[i] = true;
			numOfOut++;
		}
	}
	inWells.close();

	MaskForWells.resize(numOfWells, true);
	if (!useInWells)
		for (int i = 0; i < numOfWells; i++) {
			if (!WellsType[i]) MaskForWells[i] = false;
		}
	if (!useOutWells)
		for (int i = 0; i < numOfWells; i++) {
			if (WellsType[i]) MaskForWells[i] = false;
		}


	ifstream inLayers("ContoursMaterials");
	inLayers >> numOfLayers;

	ifstream inContours("ContoursMaterialGeometry");
	inContours >> numOfContours;
	for (int i = 0; i < numOfContours; i++) {
		vector<int> tempV;

		inContours >> tempI;
		tempV.push_back(numOfLayers - tempI);
		inContours >> tempI;
		tempV.push_back(numOfLayers - tempI);
		ContourBound.push_back(tempV);

		inContours >> temp;
		inContours >> temp;
		inContours >> temp;
		inContours >> tempI;
		for (int j = 0; j < tempI; j++) {
			inContours >> temp;
			inContours >> temp;
		}
	}
	inContours.close();

	for (int i = 0; i < numOfLayers; i++) {

		inLayers >> temp;
		inLayers >> tempD;

		TrueParams.push_back(tempD * scale); //K

		inLayers >> temp;
		inLayers >> tempD;

		TrueParams.push_back(tempD); //Phi

		for (int j = 0; j < 10; j++) {
			inLayers >> temp;
		}
	}
	for (int i = 0; i < numOfContours; i++) {

		inLayers >> temp;
		inLayers >> temp;
		inLayers >> temp;
		inLayers >> tempD;

		TrueParams.push_back(tempD * scale); //K

		inLayers >> temp;
		inLayers >> tempD;

		TrueParams.push_back(tempD); //Phi

		for (int j = 0; j < 10; j++) {
			inLayers >> temp;
		}
	}
	inLayers.close();
}

void GaussNewton::copyFiles(string from, string to) {
	if (!copyAll) {
		for (int i = 0; i < numOfWells; i++)
			if (MaskForWells[i])
				for (int j = 0; j < numOfPlots; j++)
					CopyFile((from + "\\" + plots[i][j]).c_str(), (to + "\\" + plots[i][j]).c_str(), FALSE);
	}
	else {
		for (int i = 1; i <= numOfWells; i++)
			if (MaskForWells[i - 1]) {
				CopyFile((from + "\\" + to_string(i) + "_s2.txt").c_str(), (to + "\\" + to_string(i) + "_s2.txt").c_str(), FALSE);
				CopyFile((from + "\\" + "avg_m_per_" + to_string(i) + "_s2.txt").c_str(), (to + "\\" + "avg_m_per_" + to_string(i) + "_s2.txt").c_str(), FALSE);
				CopyFile((from + "\\" + to_string(i) + "_P.txt").c_str(), (to + "\\" + to_string(i) + "_P.txt").c_str(), FALSE);
			}
	}
}

void GaussNewton::SetReverse() {

	MaskForLayers.resize(numOfLayers + numOfContours, true);

	cout << "\nModel\n - " << numOfLayers << " layers\n - " << numOfContours << " contours\n - ";
	cout << numOfWells << " wells (In: " << numOfIn << ", Out: " << numOfOut << ")\n\n";
	int x;
	offLayerSize = 0;
	offContoursSize = 0;
	vector<int> offLayers, offContours;

	if (numOfLayers + numOfContours < 2) goto next;
	cout << "Select nonsuitable layers from: ";
	for (int i = 1; i <= numOfLayers; i++) cout << i << " ";

again:
	scanf_s("%c", &x); //scanf \n after model name
	cout << "\n>> ";
	offLayers.clear();
	while (cin.peek() != '\n') {
		cin >> x;
		offLayers.push_back(numOfLayers - x + 1); //слои идут в обратном порядке
	}
	offLayerSize = offLayers.size();

	if (offLayerSize == numOfLayers && numOfContours == 0) {
		cout << "You cant choose all layers";
		goto again;
	}
	for (int i = 0; i < offLayerSize; i++) {
		if (offLayers[i] > numOfLayers)
		{
			cout << "You choosed non exist layer";
			goto again;
		}
	}
	for (int i = 0; i < offLayerSize; i++) MaskForLayers[offLayers[i] - 1] = false;

next:

	if (numOfContours < 1) goto next2;
	cout << "\nSelect nonsuitable contours from:\n";
	for (int i = 1; i <= numOfContours; i++) cout << " " << i << " (T: " << ContourBound[i - 1][1]
		<< ", B: " << ContourBound[i - 1][0] << ")\n";

again2:
	scanf_s("%c", &x); //scanf \n
	cout << "\n>> ";
	offContours.clear();

	while (cin.peek() != '\n') {
		cin >> x;
		offContours.push_back(x);
	}
	offContoursSize = offContours.size();

	if (offContoursSize == numOfContours && numOfLayers == offLayerSize) {
		cout << "You cant choose all contours";
		goto again2;
	}
	for (int i = 0; i < offContoursSize; i++) {
		if (offContours[i] > numOfContours)
		{
			cout << "You choosed non exist contour";
			goto again2;
		}
	}
	for (int i = 0; i < offContoursSize; i++) MaskForLayers[offContours[i] - 1 + numOfLayers] = false;

next2:

	N = numOfPlots * timePoints * numOfWells;
	M = 2 * (numOfLayers + numOfContours);
	n = numOfPlots * timePoints * (useInWells * numOfIn + useOutWells * numOfOut);
	m = 2 * ((numOfLayers + numOfContours) - (offLayerSize + offContoursSize));
}

void GaussNewton::readFx(Vector& fx, string pp) {
	fx.clear();
	string temp;
	_chdir((pp + "\\output").c_str());
	for (int i = 0; i < numOfWells; i++) {
		if (MaskForWells[i]) {
			for (int j = 0; j < numOfPlots; j++) {
				ifstream file(pp + "\\output\\" + plots[i][j]);
				for (int k = 0; k < timePoints; k++) {
					file >> temp;
					file >> temp;
					fx.push_back((atof(temp.c_str()) == 0) ? eps : atof(temp.c_str()));
				}
				file.close();
			}
		}
	}

}

void GaussNewton::setParameters(Vector& params, string pp) {
	_chdir((pp + "\\mesh").c_str());

	int n1 = 3, n2 = 5, nn = 1, i = 0, j = 0, maskBlock = 0;

	ifstream in(pp + "\\mesh\\ContoursMaterials");
	ofstream out(pp + "\\mesh\\ContoursMaterials_temp");

	if (!in || !out)
	{
		cout << "Error opening files!" << endl;
		return;
	}
	string strTemp;

	while (nn <= numOfLayers * 14 + 1) {
		in >> strTemp;
		if (nn == n1) {
			if (MaskForLayers[maskBlock])
			{
				out << scientific << std::fixed << std::setprecision(22) << params[j] / scale << "\n";
				j++;
			}
			else out << scientific << std::fixed << std::setprecision(22) << TrueParams[i] / scale << "\n";
			n1 += 14;
			i++;
		}
		else {
			if (nn == n2) {
				if (MaskForLayers[maskBlock])
				{
					out << scientific << std::fixed << std::setprecision(22) << params[j] << "\n";
					j++;
				}
				else out << scientific << std::fixed << std::setprecision(22) << TrueParams[i] << "\n";
				n2 += 14;
				i++;
				maskBlock++;
			}
			else {
				out << strTemp << "\n";
			}
		}
		nn++;
	}
	n1 = (numOfLayers * 14 + 1) + 4;
	n2 = (numOfLayers * 14 + 1) + 6;
	while (nn <= (numOfLayers * 14 + 1) + 16 * numOfContours) {
		in >> strTemp;
		if (nn == n1) {
			if (MaskForLayers[maskBlock])
			{
				out << scientific << std::fixed << std::setprecision(22) << params[j] / scale << "\n";
				j++;
			}
			else out << scientific << std::fixed << std::setprecision(22) << TrueParams[i] / scale << "\n";
			n1 += 16;
			i++;
		}
		else {
			if (nn == n2) {
				if (MaskForLayers[maskBlock])
				{
					out << scientific << std::fixed << std::setprecision(22) << params[j] << "\n";
					j++;
				}
				else out << scientific << std::fixed << std::setprecision(22) << TrueParams[i] << "\n";
				n2 += 16;
				i++;
				maskBlock++;
			}
			else {
				out << strTemp << "\n";
			}
		}
		nn++;
	}
	in.close();
	out.close();

	if (seeOutput) {
		if (remove("ContoursMaterials") == 0)
			cout << "Deleted successfully: " << pp << endl;
		else {
			cout << "Unable to delete the file: " << pp << endl;
		}

		if (rename("ContoursMaterials_temp", "ContoursMaterials") != 0)
			cout << "Error renaming file: " << pp << endl;
		else
			cout << "File renamed successfully: " << pp << endl;
	}
	else {
		remove("ContoursMaterials");
		rename("ContoursMaterials_temp", "ContoursMaterials");
	}
}

void GaussNewton::init(double k, double phi) {

	trash = (seeOutput) ? "" : "> trash.txt";

	_chdir(pathMesh.c_str());
	system(("cd " + pathMesh + " & C:\\Filtr\\MeshMatContours.exe" + trash).c_str()); //формирование среды

	ReadData();

	system(("cd " + path + " & C:\\Filtr\\Filtr3D.exe" + trash).c_str()); //прямой решатель

	//формируем набор названий файлов для использующихся графиков 
	numOfPlots = useS + useM + useP;
	for (int i = 1; i <= numOfWells; i++) {
		vector<string> tempS;
		if (useS) tempS.push_back(to_string(i) + "_s2.txt");
		if (useM) tempS.push_back("avg_m_per_" + to_string(i) + "_s2.txt");
		if (useP) tempS.push_back(to_string(i) + "_P.txt");
		plots.push_back(tempS);
	}

	_mkdir(pathReverse.c_str());
	string pp = pathReverse + "\\T";
	_mkdir(pp.c_str());
	copyFiles(pathOutput, pp);

	_chdir(pathOutput3D.c_str());
	ifstream ff("times_main");
	ff >> timePoints;
	ff.close();
	timePoints--; //в моментальном отборе последней точки не хватает 

	SetReverse(); //выясняем в каких слоях и контурах будет подбор

	if (PARALLEL) {
		if (m <= Threads) Threads = m;
		else if (m % Threads != 0) {
			while (m % Threads != 0) Threads--;

		}
		omp_set_num_threads(Threads);

		_chdir(pathHome.c_str());
#pragma omp parallel
		{
			int ID = omp_get_thread_num();
			string th = pathHome + "\\thread_" + to_string(ID);
			std::filesystem::copy(path, th, std::filesystem::copy_options::recursive);

#pragma omp barrier
		}
	}
	if (customFx) {
		int i;
		cout << "Put your Fx, and enter any symbol\n";
		cin >> i;
	}
	readFx(TrueFx, path);

	//p_init просто чтобы сохранить начальное приближение 
	for (int i = 0; i < m / 2; i++) {

		Params.push_back(k * scale); //K
		p_init.push_back(k * scale); //K
		Params.push_back(phi); //Phi
		p_init.push_back(phi); //Phi
	}
	setParameters(Params, path);
}

inline void GaussNewton::f(Vector& params, Vector &fx) {

	setParameters(params, path);

	_chdir(pathMesh.c_str());
	system(("cd " + pathMesh + " & C:\\Filtr\\MeshMatContours.exe" + trash).c_str()); //формирование среды
	
	
	_chdir(path.c_str());
	system(("cd " + path + " & C:\\Filtr\\Filtr3D.exe" + trash).c_str()); //прямой решатель

	readFx(fx, path);
}

double GaussNewton::Functional(Vector& truefx, Vector& fx) {
	double F = 0, f;
#pragma omp parallel for private(f) reduction(+:F)
	for (int i = 0; i < n; i++) {
		f = (truefx[i] - fx[i]) / fx[i];
		F += f * f;
	}
	return F;

}

void GaussNewton::FunctionalSeparate(Vector& truefx, Vector& fx, int iter) {
	Vector F(numOfPlots,0);
	double f;
	int flag;
	for (int i = 0; i < n; i++) {
		flag = ((int)(i / timePoints)) % numOfPlots;
		f = (truefx[i] - fx[i]) / fx[i];
		F[flag] += f * f;
	}
	for (int i = 0; i < numOfPlots; i++) {
		FunctionalProgress[iter].push_back(F[i]);

	}

}

void GaussNewton::calcJ() {
	// Производные
	double dp = 1.1;
#pragma omp parallel for
	for (int j = 0; j < m; j++)
		if (abs(p0[j]) < eps) p0[j] = sign(p0[j]) * eps;
	int DIV = m / Threads;
	for (int k = 0; k < DIV; k++)
	{
		for (int j = k * Threads; j < (k + 1) * Threads; j++) {
			p0[j] *= dp;
			string pp1 = pathHome + "\\thread_" + to_string(j - k * Threads);
			setParameters(p0, pp1);

			_chdir((pp1 + "\\mesh").c_str());
			system(("cd " + pp1 + "\\mesh & C:\\Filtr\\MeshMatContours.exe" + trash).c_str()); //формирование среды

			p0[j] /= dp;
		}
#pragma omp parallel
		{
			int ID = omp_get_thread_num();
			string pp2 = pathHome + "\\thread_" + to_string(ID);
			Sleep(100 * ID);
			system(("cd " + pp2 + " & C:\\Filtr\\Filtr3D.exe" + trash).c_str()); //прямой решатель
#pragma omp barrier
		}
		for (int j = k * Threads; j < (k + 1) * Threads; j++)
		{
			string pp3 = pathHome + "\\thread_" + to_string(j - k * Threads);
			Vector Fx1(n);
			readFx(Fx1, pp3);
#pragma omp parallel for
			for (int i = 0; i < n; ++i)
				J[j][i] = (Fx[i] - Fx1[i]) / ((dp - 1) * p0[j]);
		}
	}
	if (seeJacobian) {
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < m; j++)
				cout << left << setw(10) << J[j][i];
			cout << endl;
		}
	}
}

Vector GaussNewton::SolveGauss(Matrix& A, Vector& b)
{
	int m = A.size();
	Vector x(m);
	Matrix B(m, Vector(m, 0.0));
#pragma omp parallel for
	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < m; j++) B[i][j] = A[i][j];
		x[i] = b[i];
	}

	for (int k = 0; k < m - 1; k++)
	{
		for (int i = k + 1; i < m; i++)
		{
			double elem = -B[i][k] / B[k][k];
			for (int j = k; j < m; j++)
			{
				B[i][j] += elem * B[k][j];
			}
			x[i] += elem * x[k];
		}
	}

	for (int i = m - 1; i >= 0; i--)
	{
		double elem = 0;
		for (int j = i + 1; j < m; j++)
		{
			elem += B[i][j] * x[j];
		}
		x[i] = (x[i] - elem) / B[i][i];
	}
	return x;
}

void GaussNewton::method() {

	Fx.resize(n);
	p0.resize(m), p.resize(m), Delta_k.resize(m), f_.resize(m), alpha.resize(m), falpha.resize(m);
	A.resize(m, Vector(m)), Aalpha.resize(m, Vector(m)), J.resize(m, Vector(n));

	// Первое приближение
#pragma omp parallel for
	for (int i = 0; i < m; ++i) {
		if (abs(Params[i]) < eps)
			p0[i] = sign(Params[i]) * eps;
		else p0[i] = Params[i];
	}

	f(p0, Fx);
	Func0 = Func = Functional(TrueFx, Fx);
	Vector tempVec;
	tempVec.push_back(Func0);
	FunctionalProgress.push_back(tempVec);
	FunctionalSeparate(TrueFx, Fx, 0);

	string pp = pathReverse + "\\0";
	_mkdir(pp.c_str());
	copyFiles(pathOutput, pp);

	// Гаусс-Ньютон
	for (size_t iter = 1; iter < MaxIters; iter++) {
		std::cout << "\nITER = " << iter << endl;
		ITERATIONS = iter;

		if(PARALLEL) calcJ();
		else {
			Vector Fx1(n);
			for (size_t j = 0; j < m; j++) {
				if (abs(p0[j]) < eps) p0[j] = sign(p0[j]) * eps;
				double dp = 1.1;
				p0[j] *= dp;
				f(p0, Fx1);
				p0[j] /= dp;
				for (size_t i = 0; i < n; ++i)
					J[j][i] = (Fx[i] - Fx1[i]) / ((dp - 1) * p0[j]);
			}
		}
		

		std::cout << '\n';

		// Решение СЛАУ
		for (size_t i = 0; i < m; ++i) {
			f_[i] = 0.;
			for (size_t k = 0; k < n; ++k)
				f_[i] -= J[i][k] * (TrueFx[k] - Fx[k]) / (Fx[k] * Fx[k]);
			for (size_t j = 0; j < m; ++j) {
				A[i][j] = 0.0;
				for (size_t k = 0; k < n; ++k)
					A[i][j] += (J[i][k] * J[j][k]) / (Fx[k] * Fx[k]);
			}
		}

		// Регуляризация
		for (size_t i = 0; i < m; ++i)
			alpha[i] = AlphaMin;
		bool stop = false;
		bool stopalpha = false;
		do {
			// СЛАУ с регуляризирующим параметром
			for (size_t i = 0; i < m; ++i) {
				for (size_t j = 0; j < m; ++j)
					Aalpha[i][j] = A[i][j];
				Aalpha[i][i] += alpha[i];
				falpha[i] = f_[i] - alpha[i] * (p0[i] - Params[i]);
			}
			//Метод Гаусса
			Delta_k = SolveGauss(Aalpha, falpha);

			// Новое приближение
			for (size_t i = 0; i < m; ++i)
				p[i] = p0[i] + Delta_k[i];

			stop = true;
			for (size_t i = 0; i < m; i += 2) 
			{
				if (p[i] < IntervalK[0] || p[i] > IntervalK[1]) {
					alpha[i] *= 1.1;
					stop = false;
				}
				if (p[i+1] < IntervalPhi[0] || p[i+1] > IntervalPhi[1]) {
					alpha[i + 1] *= 1.1;
					stop = false;
				}
			}

			double alphaMax = 0.;
			for (size_t i = 0; i < m; ++i)
				if (alpha[i] > alphaMax)
					alphaMax = alpha[i];
			for (size_t i = 0; i < m; ++i)
				while (alphaMax / alpha[i] > 1e5) {
					alpha[i] *= 1.1;
				}
			// Проверка на размерность
			for (size_t i = 0; i < m && !stopalpha; ++i)
				if (alpha[i] > AlphaMax) {
					std::cout << "Exit by max alpha" << std::endl;
					stopalpha = true;
					stop = true;
				}
		} while (!stop);

		std::cout << "Alpha = ";
		for (size_t i = 0; i < m; ++i)
			std::cout << alpha[i] << ' ';
		std::cout << '\n';

		std::cout << "Parameters = ";
		for (size_t i = 0; i < m; ++i)
			std::cout << p[i] << ' ';
		std::cout << '\n';

		// Релаксация
		double beta = 1.;
		f(p, Fx);

		double Func1 = Functional(TrueFx, Fx);
		stop = false;
		while (Func1 >= Func) {
			std::cout << "  beta = " << beta << '\n';
			beta /= 2;
#pragma omp parallel for
			for (int i = 0; i < m; ++i)
				p[i] = p0[i] + beta * Delta_k[i];
			f(p, Fx);
			Func1 = Functional(TrueFx, Fx);
			if (beta < BetaMin) {
				stop = true;
				break;
			}
		}

		// Обновление параметров
		Func = Func1;
		for (size_t i = 0; i < m; ++i) {
			p0[i] = p[i];
			Params[i] = p0[i];
		}
		std::cout << "Func/Func0 = " << Func / Func0 << '\n';

		string pp = pathReverse + "\\" + to_string(ITERATIONS);
		_mkdir(pp.c_str());
		copyFiles(pathOutput, pp);

		Vector tempVec;
		tempVec.push_back(Func);
		FunctionalProgress.push_back(tempVec);
		FunctionalSeparate(TrueFx, Fx, ITERATIONS);

		// Условие выхода
		if (Func / Func0 < Eps ||  stop) break;
		for (int i = 0; i < m; i++)
			std::cout << p0[i] << ' ';
		std::cout << '\n';

	}
	if (PARALLEL) {
		#pragma omp parallel
		{
			int ID = omp_get_thread_num();
			string th = pathHome + "\\thread_" + to_string(ID);
			std::filesystem::remove_all(th);
			
		#pragma omp barrier
		}
	}
	ReturnTrue(TrueParams, path); //восстановление файла материалов для следующего запуска
}

void GaussNewton::Results() {

	Vector new_p;
	for (int i = 2 * (numOfLayers - offLayerSize) - 1; i > 0; i = i - 2) {
		new_p.push_back(p[i - 1]);
		new_p.push_back(p[i]);
	}
	for (int i = 2*(numOfLayers - offLayerSize); i < p.size(); i = i + 2) {
		new_p.push_back(p[i]);
		new_p.push_back(p[i+1]);
	}
	Vector new_TrueParams;
	for (int i = 2 * numOfLayers - 1; i > 0; i = i - 2) {
		new_TrueParams.push_back(TrueParams[i - 1]);
		new_TrueParams.push_back(TrueParams[i]);
	}
	for (int i = 2 * numOfLayers; i < TrueParams.size(); i = i + 2) {
		new_TrueParams.push_back(TrueParams[i]);
		new_TrueParams.push_back(TrueParams[i + 1]);
	}
	if (p.size() != new_p.size() || TrueParams.size() != new_TrueParams.size()) {
		cout << "ERROR reversing vector\n";
		return;
	}
	ofstream results(pathReverse + "\\results.txt");
	int flag = 0;
	results << "\n\n\tStart\tResult\tTrue\n\nLayers\n";
	cout << "\n\n" << left << setw(9) << " " << setw(12) << "Start" << setw(12) << "Result" << setw(12) << "True\n\nLayers\n";
	for (int i = 0, j = 0, k; i < M; i++) {
		k = (i % 2 == 0) ? 1E15 / scale : 1;
		if (i == 2 * numOfLayers) {
			results << "\nContours\n";
			cout << "\nContours" << endl;
			flag = 1; }
		if (MaskForLayers[(int)i / 2])
		{
			results << ((i % 2 == 0) ? to_string(i / 2 - flag * numOfLayers + 1) : "") << "\t" << to_string(p_init[j] * k) << "\t" << to_string(new_p[j] * k)
				<< "\t" << to_string(new_TrueParams[i] * k) << "\n";
			cout << left << setw(9) << ((i % 2 == 0) ? to_string(i/2 - flag * numOfLayers + 1) : "")
				<< setw(12) << p_init[j] * k 
				<< setw(12) << new_p[j] * k 
				<< setw(12) << new_TrueParams[i] * k << endl;
			j++;
		}
		else {
			results << ((i % 2 == 0) ? to_string(i / 2 - flag * numOfLayers + 1) : "")
				<< "\t" << to_string(new_TrueParams[i] * k)
				<< "\t" << to_string(new_TrueParams[i] * k)
				<< "\t" << to_string(new_TrueParams[i] * k) << "\n";
			cout << left << setw(9) << ((i % 2 == 0) ? to_string(i / 2 - flag * numOfLayers + 1) : "")
				<< setw(12) << new_TrueParams[i] * k
				<< setw(12) << new_TrueParams[i] * k
				<< setw(12) << new_TrueParams[i] * k << endl;
		}
	}
	results << "\nIterations: " << ITERATIONS << "\n";
	results << "Functional: " << Func << "\n";
	cout << "\nIterations: " << ITERATIONS << endl;
	cout << "Functional: " << Func << endl;

	results.close();

	ofstream results2(pathReverse + "\\FunctionalProgress.txt");
	results2 << "Iter\tF\t";
	if (useS) results2 << "F(S)\t";
	if (useM) results2 << "F(M)\t";
	if (useP) results2 << "F(P)\t";
	for (int i = 0; i < FunctionalProgress.size(); i++) 
	{
		results2 << "\n" << to_string(i) << "\t";
		for (int j = 0; j < FunctionalProgress[0].size(); j++) 		
			results2 << FunctionalProgress[i][j] << "\t";		
		}
	results2.close();
}

void GaussNewton::ReturnTrue(Vector& params, string pp) {
	_chdir((pp + "\\mesh").c_str());

	int n1 = 3, n2 = 5, nn = 1, i = 0;

	ifstream in(pp + "\\mesh\\ContoursMaterials");
	ofstream out(pp + "\\mesh\\ContoursMaterials_temp");

	if (!in || !out)
	{
		cout << "Error opening files!" << endl;
		return;
	}
	string strTemp;

	while (nn <= numOfLayers * 14 + 1) {
		in >> strTemp;
		if (nn == n1) {
			 out << scientific << std::fixed << std::setprecision(22) << TrueParams[i] / scale << "\n";
			n1 += 14;
			i++;
		}
		else {
			if (nn == n2) {
				out << scientific << std::fixed << std::setprecision(22) << TrueParams[i] << "\n";
				n2 += 14;
				i++;
			}
			else {
				out << strTemp << "\n";
			}
		}
		nn++;
	}
	n1 = (numOfLayers * 14 + 1) + 4;
	n2 = (numOfLayers * 14 + 1) + 6;
	while (nn <= (numOfLayers * 14 + 1) + 16 * numOfContours) {
		in >> strTemp;
		if (nn == n1) {
			out << scientific << std::fixed << std::setprecision(22) << TrueParams[i] / scale << "\n";
			n1 += 16;
			i++;
		}
		else {
			if (nn == n2) {
				out << scientific << std::fixed << std::setprecision(22) << TrueParams[i] << "\n";
				n2 += 16;
				i++;
			}
			else {
				out << strTemp << "\n";
			}
		}
		nn++;
	}
	in.close();
	out.close();


	if (seeOutput) {
		if (remove("ContoursMaterials") == 0)
			cout << "Deleted successfully: " << pp << endl;
		else {
			cout << "Unable to delete the file: " << pp << endl;
		}

		if (rename("ContoursMaterials_temp", "ContoursMaterials") != 0)
			cout << "Error renaming file: " << pp << endl;
		else
			cout << "File renamed successfully: " << pp << endl;
	}
	else {
		remove("ContoursMaterials");
		rename("ContoursMaterials_temp", "ContoursMaterials");
	}
}
