#include "mpi.h"
#include <stdio.h>
#include <malloc.h>
#include <math.h>
#include <time.h>
#include <stdlib.h>

// Định nghĩa các tag trao đổi dữ liệu
#define TAG_ABOVE_BOUNDARY 1
#define TAG_BELOW_BOUNDARY 2
#define TAG_STOP 3
#define TAG_CONCENTRATION 4
#define TAG_BACILLUS 5
#define TAG_WORKLOAD 6
#define TAG_ERROR 7
#define TAG_SUM_CONCETRATION 8

// 2 mảng hằng thể hiện mức chênh lệch hàng và chênh lệch cột giữa 1 ô với các ô kề cạnh của nó
const int dl[] = {1, 0, -1, 0};
const int dm[] = {0, 1, 0, -1};

// Kích thước lưới
int N = 80;
// Số tế bào trực khuẩn ban đầu
int numInitialBacillus;
// Danh sách các ô xuất phát của trực khuẩn
struct InitialBacillus {
	int l;
	int m;
} *initialBacillus;
// Hệ số phát triển
double eta = 1.00;
// Hệ số mixing
double omega = 1.50;
// Điều kiện dừng thủ tục lặp để giải phương trình khuếch tán
double epsilon = 1e-3;
// Số vô cùng nhỏ. Mọi số nhỏ hơn số này được coi như bằng 0
double infsm = 1e-6;
// Số vòng lặp phát triển của trực khuẩn
int numIteration = 1000;

// Đọc các thông số cấu hình mô phỏng
void readProblemConfiguration() {
	FILE * f = fopen("config.txt", "r");

	fscanf(f, "%d", &N);

	// Lấy thông tin về trạng thái khởi tạo của trực khuẩn
	fscanf(f, "%d", &numInitialBacillus);
	int i;
	initialBacillus = (struct InitialBacillus *) malloc (numInitialBacillus * sizeof (struct InitialBacillus));
	for (i = 0; i < numInitialBacillus; i++) {
		fscanf(f, "%d%d", &(initialBacillus[i].l), &(initialBacillus[i].m));
	}

	fscanf(f, "%lf", &eta);
	eta = round(eta * 1e2) / (1e2);

	fscanf(f, "%lf", &omega);
	omega = round(omega * 1e2) / (1e2);

	fscanf(f, "%lf", &epsilon);
	epsilon = round(epsilon * 1e5) / (1e5);

	fscanf(f, "%lf", &infsm);
	infsm = round(infsm * 1e8) / (1e8);

	fscanf(f, "%d", &numIteration);

	fclose(f);
}

// Tính giá trị lớn hơn giữa 2 số double
double max(double a, double b) {
	return a > b ? a : b;
}

// Tính giá trị lớn hơn giữa 2 số int
double maxInt(int a, int b) {
	return a > b ? a : b;
}

// Công thức cập nhật SOR
double _SOR(double *concentration, int l, int m) {
	return fabs((omega / 4) * (
		concentration[N * (l - 1) + m] +
		concentration[N * (l + 1) + m] +
		concentration[N * l + ((m + 1 + N) % N)] +
		concentration[N * l + ((m - 1 + N) % N)]) +
		(1 - omega) * concentration[N * l + m]);
}

// Hàm log lưu lại trạng thái lưới sau mỗi vòng lặp phát triển
void iteratorLog(int iter, double *concentration, int *bacillus) {
	char *logFile = (char *) malloc (1024 * sizeof(char));
	sprintf(logFile, "./log/log_N=%d_eta=%0.2lf_omega=%0.2lf_iter=%d.log", N, eta, omega, iter);
	FILE *f = fopen(logFile, "w");
	int l, m;
	for (l = 0; l < N; l++) {
		for (m = 0; m < N; m++) {
			fprintf(f, "%lf ", concentration[N * l + m]);
		}
		fprintf(f, "\n");
	}
	fprintf(f, "\n");
	for (l = 0; l < N; l++) {
		for (m = 0; m < N; m++) {
			fprintf(f, "%d ", bacillus[N * l + m]);
		}
		fprintf(f, "\n");
	}
	printf("%s\n", logFile);
	fclose(f);
	free(logFile);
}

// Các process rank > 0 trao đổi thông tin về lượng thức ăn tại biên của mình
void boundariesConcentrationExchange(int rank, int size,
	double *concentration, int *bacillus, int *workload) {
	
	MPI_Status status;
	int rWorkload = workload[rank] - workload[rank - 1];

	MPI_Request request[4];

	if (rank > 1) {
		MPI_Isend(concentration + N, N, MPI_DOUBLE,
			rank - 1, TAG_ABOVE_BOUNDARY, MPI_COMM_WORLD, &(request[0]));
	}
	if (rank < size - 1) {
		MPI_Isend(concentration + N * rWorkload, N, MPI_DOUBLE,
			rank + 1, TAG_BELOW_BOUNDARY, MPI_COMM_WORLD, &(request[1]));
	}
	if (rank < size - 1) {
		MPI_Irecv(concentration + N * (rWorkload + 1), N, MPI_DOUBLE,
			rank + 1, TAG_ABOVE_BOUNDARY, MPI_COMM_WORLD, &(request[2]));
	}
	if (rank > 1) {
		MPI_Irecv(concentration, N, MPI_DOUBLE,
			rank - 1, TAG_BELOW_BOUNDARY, MPI_COMM_WORLD, &(request[3]));
	}

	if (rank < size - 1) {
		MPI_Wait(&(request[2]), &status);
	}
	if (rank > 1) {
		MPI_Wait(&(request[3]), &status);
	}
}

// Các process rank > 0 trao đổi thông tin về trực khuẩn tại biên của mình
void boundariesBacillusExchange(int rank, int size,
	double *concentration, int *bacillus, int *workload) {

	int rWorkload = workload[rank] - workload[rank - 1];
	MPI_Status status;

	MPI_Request request[4];

	if (rank > 1) {
		MPI_Isend(bacillus + N, N, MPI_INT,
			rank - 1, TAG_ABOVE_BOUNDARY, MPI_COMM_WORLD, &(request[0]));
	}
	if (rank < size - 1) {
		MPI_Isend(bacillus + N * rWorkload, N, MPI_INT,
			rank + 1, TAG_BELOW_BOUNDARY, MPI_COMM_WORLD, &(request[1]));
	}
	if (rank < size - 1) {
		MPI_Irecv(bacillus + N * (rWorkload + 1), N, MPI_INT,
			rank + 1, TAG_ABOVE_BOUNDARY, MPI_COMM_WORLD, &(request[2]));
	}
	if (rank > 1) {
		MPI_Irecv(bacillus, N, MPI_INT,
			rank - 1, TAG_BELOW_BOUNDARY, MPI_COMM_WORLD, &(request[3]));
	}

	if (rank < size - 1) {
		MPI_Wait(&(request[2]), &status);
	}
	if (rank > 1) {
		MPI_Wait(&(request[3]), &status);
	}
}

// Cập nhật mật độ thức ăn
// Thao tác này là 1 lần lặp trong thuật toán giải phương trình khuếch tán
// rb = 0 tương ứng với thao tác cập nhật ô đỏ
// rb = 1 tương ứng với thao tác cập nhật ô đen
double updateConcentration(int rank, int size,
	double *concentration, int *bacillus, int *workload, int rb) {

	double maxError = 0;
	boundariesConcentrationExchange(rank, size, concentration, bacillus, workload);

	int l, m;
	int rWorkload = workload[rank] - workload[rank - 1];
	for (l = 1; l <= rWorkload; l++) {
		for (m = 0; m < N; m++) {
			// color: màu sắc của ô đang xét
			int color = (l - 1) + m + N * workload[rank - 1];
			if (color % 2 == rb) {
				double oldConcentration = concentration[N * l + m];
				concentration[N * l + m] = _SOR(concentration, l, m);
				if (bacillus[N * l + m]) {
					concentration[N * l + m] = 0;
				}
				if (rank == 1 && l == 1) {
					concentration[N * l + m] = 0;
				}
				if (rank == size - 1 && l == rWorkload) {
					concentration[N * l + m] = 1;
				}
				maxError = max(maxError, fabs(concentration[N * l + m] - oldConcentration));
			}
		}
	}
	return maxError;
}

// Process rank 0 khởi tạo và phân chia công việc
void initialize(int rank, int size, double **concentration, int **bacillus, int **workload) {

	readProblemConfiguration();

	// Mảng concentration chứa thông tin về mật độ thức ăn trên toàn lưới
	*concentration = (double *) malloc (N * N * sizeof(double));
	// Mảng bacillus chứa thông tin về trực khuẩn trên toàn lưới
	*bacillus = (int *) malloc (N * N * sizeof(int));
	// Mảng workload chứa thông tin về khối lượng công việc mà mỗi process phải hoàn thành
	// workoad[i] = tổng số hàng mà các process rank <= i phải thực hiện
	*workload = (int *) malloc (size * sizeof(int));

	// Khởi tạo trạng thái lưới: Trừ source có c = 1, tất cả các ô khác đều có c = 0
	int l, m;
	for (m = 0; m < N; m++) {
		for (l = 0; l < N - 1; l++) {
			(*concentration)[N * l + m] = 0;
			(*bacillus)[N * l + m] = 0;
		}
		(*concentration)[N * (N - 1) + m] = 1;
		(*bacillus)[N * (N - 1) + m] = 0;
	}
	int i;
	// Khởi tạo trạng thái trực khuẩn theo mảng initialBacillus đọc được từ file config.txt
	for (i = 0; i < numInitialBacillus; i++) {
		l = initialBacillus[i].l;
		m = initialBacillus[i].m;
		(*bacillus)[N * l + m] = 1;
	}

	// Phân chia công việc cho các process rank > 0
	// avg_workload là lượng công việc ít nhất mà mỗi process rank > 0 phải thực hiện
	int avg_workload = N / (size - 1);
	// over_workload là số process phải tính toán thêm 1 hàng so với mức tối thiểu
	int over_workload = N % (size - 1);

	(*workload)[0] = 0;
	int r;
	for (r = 1; r < size; r++) {
		(*workload)[r] = (*workload)[r - 1] + avg_workload;
		if (r <= over_workload) {
			(*workload)[r] += 1;
		}
	}
}

// Process rank 0 gửi thông tin công việc đến các process khác
void divideTask(int rank, int size, double *concentration, int *bacillus, int *workload) {
	int r;
	MPI_Request request[4];
	for (r = 1; r < size; r++) {
		int rWorkload = workload[r] - workload[r - 1];
		MPI_Isend(workload, size, MPI_INT,
			r, TAG_WORKLOAD, MPI_COMM_WORLD, &(request[0]));
		MPI_Isend(concentration + N * workload[r - 1], N * rWorkload, MPI_DOUBLE,
			r, TAG_CONCENTRATION, MPI_COMM_WORLD, &(request[1]));
		MPI_Isend(bacillus + N * workload[r - 1], N * rWorkload, MPI_INT,
			r, TAG_BACILLUS, MPI_COMM_WORLD, &(request[2]));
	}
}

// Các process rank > 0 nhận việc từ process rank 0
void recvTask(int rank, int size, double **concentration, int **bacillus, int **workload) {

	MPI_Status status;
	MPI_Request request[3];

	*workload = (int *) malloc (size * sizeof(int));
	MPI_Irecv(*workload, size, MPI_INT,
		0, TAG_WORKLOAD, MPI_COMM_WORLD, &(request[0]));
	MPI_Wait(&(request[0]), &status);

	// rWorkdload là lượng công việc mà process hiện tại phải thực hiện
	int rWorkload = (*workload)[rank] - (*workload)[rank - 1];
	
	// concentration là lượng thức ăn trên các hàng mà process hiện tại quản lý
	// concentration đồng thời chứa thông tin các hàng biên của 2 process lân cận
	*concentration = (double *) malloc (N * (rWorkload + 2) * sizeof(double));
	MPI_Irecv((*concentration) + N, N * rWorkload, MPI_DOUBLE,
		0, TAG_CONCENTRATION, MPI_COMM_WORLD, &(request[1]));

	// bacillus là thông tin về trực khuẩn trên các hàng mà process hiện tại quản lý
	// bacillus đồng thời chứa thông tin các hàng biên của 2 process lân cận
	*bacillus = (int *) malloc (N * (rWorkload + 2) * sizeof(int));
	MPI_Irecv((*bacillus) + N, N * rWorkload, MPI_INT,
		0, TAG_BACILLUS, MPI_COMM_WORLD, &(request[2]));

	if (rank == 1) {
		int m;
		for (m = 0; m < N; m++) {
			(*bacillus)[m] = 0;
		}
	} else if (rank == size - 1) {
		int m;
		for (m = 0; m < N; m++) {
			(*bacillus)[N * (rWorkload + 1) + m] = 0;
		}
	}

	MPI_Wait(&(request[1]), &status);
	MPI_Wait(&(request[2]), &status);
}

// Process rank 0 thu thập thông tin về thức ăn trên toàn lưới và kiểm tra điều kiện dừng
void collectConcentration(int rank, int size,
	double *concentration, int *bacillus, int *workload) {

	// Cực đại lỗi cục bộ
	double error;
	// Cực đại lỗi toàn cục
	double maxError;
	// Cờ kết thúc cho biết điều kiện dừng đã thỏa mãn hay chưa
	int stop = 0;
	int i = 0;
	MPI_Status status;
	MPI_Request request[4];

	int r;
	while (stop != 1) {
		for (r = 1; r < size; r++) {
			MPI_Isend(&stop, 1, MPI_INT,
				r, TAG_STOP, MPI_COMM_WORLD, &(request[2]));
		}
		maxError = 0;
		i++;
		for (r = 1; r < size; r++) {
			int rWorkload = workload[r] - workload[r - 1];
			MPI_Irecv(&error, 1, MPI_DOUBLE,
				r, TAG_ERROR, MPI_COMM_WORLD, &(request[1]));
			MPI_Wait(&(request[1]), &status);
			maxError = max(maxError, error);
		}
		// Khi cực đại lỗi toàn cục < epsilon, điều kiện dừng được thỏa mãn
		if (maxError < epsilon) {
			stop = 1;
		}
	}
	// Gửi tín hiệu điều khiển vòng lặp (dừng hay tiếp tục) cho các process rank > 0
	// Nhận các dữ liệu cuối cùng của các process rank > 0
	for (r = 1; r < size; r++) {
		MPI_Isend(&stop, 1, MPI_INT,
			r, TAG_STOP, MPI_COMM_WORLD, &(request[2]));
		int rWorkload = workload[r] - workload[r - 1];
		MPI_Irecv(concentration + N * workload[r - 1], N * rWorkload, MPI_DOUBLE,
			r, TAG_CONCENTRATION, MPI_COMM_WORLD, &(request[0]));
		MPI_Irecv(bacillus + N * workload[r - 1], N * rWorkload, MPI_INT,
			r, TAG_BACILLUS, MPI_COMM_WORLD, &(request[3]));
		MPI_Wait(&(request[0]), &status);
		MPI_Wait(&(request[3]), &status);
	}
}

// Giải phương trình khuếch tán
void solveDiffusionEquation(int rank, int size,
	double *concentration, int *bacillus, int *workload) {

	// Cờ tín hiệu kết thúc
	int stop = 0;
	int i = 0;
	MPI_Status status;
	MPI_Request request[4];
	// Khối lượng công việc process hiện tại cần thực hiện
	int rWorkload = workload[rank] - workload[rank - 1];

	while (stop != 1) {
		MPI_Irecv(&stop, 1, MPI_INT,
			0, TAG_STOP, MPI_COMM_WORLD, &(request[0]));
		MPI_Wait(&(request[0]), &status);
		if (stop == 1) {
			MPI_Isend(concentration + N, N * rWorkload, MPI_DOUBLE,
				0, TAG_CONCENTRATION, MPI_COMM_WORLD, &(request[1]));
			MPI_Isend(bacillus + N, N * rWorkload, MPI_INT,
				0, TAG_BACILLUS, MPI_COMM_WORLD, &(request[3]));
			break;
		}
		i++;

		// Cập nhật ô đỏ
		double maxRedError = updateConcentration(rank, size,
			concentration, bacillus, workload, 0);
		// Cập nhật ô đen
		double maxBlackError = updateConcentration(rank, size,
			concentration, bacillus, workload, 1);

		// Tính toán và gửi cực đại lỗi cục bộ cho process rank 0
		double maxError = max(maxRedError, maxBlackError);

		MPI_Isend(&maxError, 1, MPI_DOUBLE,
			0, TAG_ERROR, MPI_COMM_WORLD, &(request[2]));
	}
}

// Tính toán tổng lượng thức ăn trên các ô lân cận với trực khuẩn
void calculateGlobalCandidateConcentration(int rank, int size,
	double *concentration, int *bacillus, int *workload) {

	// Tổng lượng thức ăn trên các ô lân cận với trực khuẩn mà mỗi process tính được
	double localCandidateConcentration;
	// Tổng lượng thức ăn trên các ô lân cận với trực khuẩn trên toàn lưới
	double globalCandidateConcentration = 0;

	MPI_Status status;
	MPI_Request request[2];

	int r;
	for (r = 1; r < size; r++) {
		MPI_Irecv(&localCandidateConcentration, 1, MPI_DOUBLE,
			r, TAG_SUM_CONCETRATION, MPI_COMM_WORLD, &(request[0]));
		MPI_Wait(&(request[0]), &status);
		globalCandidateConcentration += localCandidateConcentration;
	}
	
	for (r = 1; r < size; r++) {
		MPI_Isend(&globalCandidateConcentration, 1, MPI_DOUBLE,
			r, TAG_SUM_CONCETRATION, MPI_COMM_WORLD, &(request[1]));
	}
}

// Hàm phát triển. Cho trực khuẩn phát triển ngẫu nhiên đến các ô lân cận, tuân theo công thức xác suất phát triển.
void randomGrow(int rank, int size,
	double *concentration, int *bacillus, int *workload) {

	boundariesBacillusExchange(rank, size,
		concentration, bacillus, workload);

	double localCandidateConcentration = 0;
	double globalCandidateConcentration;

	int rWorkload = workload[rank] - workload[rank - 1];
	MPI_Status status;
	MPI_Request request[2];

	int l, m, i;
	// isCand: xác định ô hiện tại có lân cận với trực khuẩn hay không
	int isCand;
	for (l = 1; l <= rWorkload; l++) {
		for (m = 0; m < N; m++) {
			isCand = 0;
			for (i = 0; i < 4; i++) {
				isCand = maxInt(isCand, bacillus[N * ((l + dl[i]) % N) + ((m + dm[i]) % N)]);
			}
			if (bacillus[N * l + m] == 0) {
				localCandidateConcentration += pow(concentration[N * l + m] * isCand, eta);
				bacillus[N * l + m] = -isCand;
			}
		}
	}

	MPI_Isend(&localCandidateConcentration, 1, MPI_DOUBLE,
		0, TAG_SUM_CONCETRATION, MPI_COMM_WORLD, &(request[0]));
	MPI_Irecv(&globalCandidateConcentration, 1, MPI_DOUBLE,
		0, TAG_SUM_CONCETRATION, MPI_COMM_WORLD, &(request[1]));
	MPI_Wait(&(request[1]), &status);

	double growProbability;
	double grow;
	for (l = 1; l <= rWorkload; l++) {
		for (m = 0; m < N; m++) {
			if (bacillus[N * l + m] == -1) {
				bacillus[N * l + m] = 0;
				// Nếu tổng lượng thức ăn trên toàn lưới quá nhỏ, ta coi như không còn thức ăn
				// Trực khuẩn k thể phát triển nếu không còn thức ăn
				if (globalCandidateConcentration <= infsm) {
					continue;
				}
				// Tính xác suất phát triển của trực khuẩn và cho trực khuẩn phát triển
				growProbability = pow(concentration[N * l + m], eta) / globalCandidateConcentration;
				grow = 1.0 * (rand() % (N * N)) / (N * N);
				if (grow < growProbability) {
					bacillus[N * l + m] = 1;
					concentration[N * l + m] = 0;
				}
			}
		}
	}
}

// Phát triển trực khuẩn
void bacillusGrow(int rank, int size, double *concentration, int *bacillus, int *workload) {

	int rWorkload = workload[rank] - workload[rank - 1];

	if (rank == 0) {
		calculateGlobalCandidateConcentration(rank, size, concentration, bacillus, workload);
	} else {
		randomGrow(rank, size, concentration, bacillus, workload);
	}
}

// Giải phương trình khuếch tán để đạt trạng thái steady
void diffusionSteady(int rank, int size,
	double *concentration, int *bacillus, int *workload, int iter) {

	if (rank == 0) {
		collectConcentration(rank, size, concentration, bacillus, workload);
	} else {
		solveDiffusionEquation(rank, size, concentration, bacillus, workload);
	}

	if (rank == 0) {
		iteratorLog(iter, concentration, bacillus);
	}

}

void main(int argc, char **argv) {

	MPI_Init(&argc, &argv);

	srand(time(0));

	int rank, size;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	double *concentration;
	int *bacillus;
	int *workload;

	if (rank == 0) {
		initialize(rank, size, &concentration, &bacillus, &workload);
		divideTask(rank, size, concentration, bacillus, workload);
	} else {
		recvTask(rank, size, &concentration, &bacillus, &workload);
	}

	int iter = 0;

	// Trạng thái steady ban đầu
	diffusionSteady(rank, size, concentration, bacillus, workload, iter);

	for (iter = 1; iter <= numIteration; iter++) {
		bacillusGrow(rank, size, concentration, bacillus, workload);
		diffusionSteady(rank, size, concentration, bacillus, workload, iter);
	}

	free(concentration);
	free(bacillus);
	free(workload);

	MPI_Finalize();

	return;
}
