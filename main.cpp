#include "reverse.h"

void main() {

	GaussNewton G;

	G.PARALLEL = true; //распараллеливание
	G.newWay = true; // true - выбор папки из консоли / false - ввод имени папки 
	G.seeOutput = false; //видеть вывод прямого решателя и работы с файлами
	G.seeJacobian = false; //видеть вывод Якобиана
	G.copyAll = true; //копировать даже неиспользуемые графики
	G.customFx = false; //подставить свои графики

	G.useOutWells = true; //использовать добывающие скважины
	G.useInWells = false; //использовать нагнетающие скважины
	G.useS = true; //использовать суммарный отбор
	G.useM = true; //использовать моментальный отбор
	G.useP = true; //использовать давление


	if (!(G.useOutWells + G.useInWells)) { cout << "Error: No well is used\n"; goto error; }
	if (!(G.useS + G.useM + G.useP)) { cout << "Error: No graphics is used\n"; goto error; }
	if (G.PARALLEL) {
		#pragma omp parallel
		{
			if (omp_get_thread_num() == 0) G.Threads = omp_get_num_threads();
		}
	}
	cout << "Threads: " << G.Threads << "\n";
	cout << (G.useS ? "S" : " ") << (G.useM ? "M" : "") << (G.useP ? "P" : "") << (G.useInWells ? " + P" : "") << "\n\n";


	G.createPath();
	G.init(3000 * 1E-15, 0.3);
	G.method();
	G.Results();

error:;
} 
