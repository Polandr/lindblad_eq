// Matrix input and output operations realization

#include <iomanip>
#include <fcntl.h>
#include <unistd.h>

using namespace std;

// Information about submatrices:

int* Matrix::gather_info()
{
	int* total_info = NULL;
	int my_info[4];
	if (ProcessorGrid::is_root())
		total_info = (int*)malloc(info.proc_grid_size()*4*sizeof(int));

	my_info[0] = info.row_offset();
	my_info[1] = info.col_offset();
	my_info[2] = n_rows;
	my_info[3] = n_cols;
	MPI_Gather(my_info,4,MPI_INT,total_info,4,MPI_INT,ProcessorGrid::root,MPI_COMM_WORLD);

	return total_info;
}

int get_proc_by_index(int* info, int row, int col)
{
	for (int i = 0; i < ProcessorGrid::size(); i++)
	{
		if (row >= info[4*i] &&
			col >= info[4*i+1] &&
			row < (info[4*i] + info[4*i+2]) &&
			col < (info[4*i+1] + info[4*i+3]))
			return i;
	}
}


// Global stream I/O------------------------------------------------------------------------

void Matrix::operator >> (std::ostream& out)
{
	out << flush;
	ProcessorGrid::barrier();
	for (int i = 0; i < global_n_rows(); i++)
		for (int j = 0; j < global_n_cols(); j++)
		{
			if (in_block(i,j))
			{
				out << setprecision(2) << setw(6) << data[(i - info.row_offset()) + n_rows * (j - info.col_offset())];
				if (j == (global_n_cols() - 1))
					out << endl;
				else
					out << ' ';
				out << flush;
			}
			ProcessorGrid::barrier();
		}
}

void Matrix::operator << (std::istream& in)
{
	ProcessorGrid::barrier();
	for (int i = 0; i < global_n_rows(); i++)
		for (int j = 0; j < global_n_cols(); j++)
		{
			if (in_block(i,j))
			{
				in >> data[(i - info.row_offset()) + n_rows * (j - info.col_offset())];
			}
			ProcessorGrid::barrier();
		}
}

void write_elems(ostream& out, double* buf, int count)
{
	for (int i = 0; i < count; i++)
	{
		out << /*setprecision(2) << */
		'(' << buf[2*i] << ',' << buf[2*i+1] << ')';
		if (i != count-1)
			out << '\t';
	}
}

void read_elems(istream& in, double* buf, int count)
{
	complexd val;
	for (int i = 0; i < count; i++)
	{
		in >> val;
		buf[2*i] = val.real();
		buf[2*i+1] = val.imag();
	}
}

void Matrix::stream_output (ostream& out)
{
	int* local_proc_info = gather_info();

	if (ProcessorGrid::is_root())
	// Root
	{
		int root_row = 0;

		for (int ofs = 0; ofs < global_n_rows()*global_n_cols(); )
		{
			int row = ofs / global_n_cols();
			int col = ofs % global_n_cols();
			int src_proc = get_proc_by_index(local_proc_info,row,col);
			int length = local_proc_info[src_proc*4+3];
			double* buf = (double*)malloc(length*2*sizeof(double));

			if (src_proc == ProcessorGrid::root)
			{
				get_row(buf,root_row);
				root_row++;
			}
			else
			{
				MPI_Status status;
				MPI_Recv(buf,2*length,MPI_DOUBLE,src_proc,0,MPI_COMM_WORLD,&status);
			}
			ofs += length;
			
			write_elems(out,buf,length);
			if (ofs % global_n_cols() == 0)
				out << endl;
			else
				out << '\t';

			free(buf);
		}

		out << flush;
	}
	else
	// Not root
	{
		double* buf = (double*)malloc(n_cols*2*sizeof(double));
		for (int i = 0; i < n_rows; i++)
		{
			get_row(buf,i);
			MPI_Send(buf,n_cols*2,MPI_DOUBLE,ProcessorGrid::root,0,MPI_COMM_WORLD);
		}
		free(buf);
	}
}

void Matrix::stream_input (istream& in)
{
	int* local_proc_info = gather_info();

	if (ProcessorGrid::is_root())
	// Root
	{
		int root_row = 0;

		for (int ofs = 0; ofs < global_n_rows()*global_n_cols(); )
		{
			int row = ofs / global_n_cols();
			int col = ofs % global_n_cols();
			int trg_proc = get_proc_by_index(local_proc_info,row,col);
			// Get number of cols of local matrix on trg_proc:
			int length = local_proc_info[trg_proc*4+3];
			double* buf = (double*)malloc(length*2*sizeof(double));

			read_elems(in,buf,length);
			ofs += length;

			if (trg_proc == ProcessorGrid::root)
			{
				set_row(buf,root_row);
				root_row++;
			}
			else
			{
				MPI_Send(buf,2*length,MPI_DOUBLE,trg_proc,0,MPI_COMM_WORLD);
			}

			free(buf);
		}
	}
	else
	// Not root
	{
		double* buf = (double*)malloc(n_cols*2*sizeof(double));
		MPI_Status status;
		for (int i = 0; i < n_rows; i++)
		{
			MPI_Recv(buf,n_cols*2,MPI_DOUBLE,ProcessorGrid::root,0,MPI_COMM_WORLD,&status);
			set_row(buf,i);
		}
		free(buf);
	}
}

ostream &operator << (ostream &out, Matrix &mat)
{
	//mat >> out;
	mat.stream_output(out);
	return out;
}

istream &operator >> (istream &in, Matrix &mat)
{
	//mat << in;
	mat.stream_input(in);
	return in;
}

// Service functions---------------------------------------------------------------------------

int bin_or_txt (const char* filename) /* 0 - error, 1 - binary, 2 - text */
{
	int fd = open(filename, O_RDONLY, 0666);
	char a[7];
	read(fd, &a, 6*sizeof(char));
	close(fd);
	a[6] = '\0';
	if (!strcmp(a, "MATRIX"))
	{	
		return 1;
	}
	else 
		if (!strcmp(a, "Matrix"))
		{
			return 2;
		}
		else 
		{
			return 0;
		}
}

// Read:

void read_header(int fd, int* dims)
{
	lseek(fd,6*sizeof(char),SEEK_CUR);
	if (read(fd, dims, 2*sizeof(int)) < 2*sizeof(int))
		throw Matrix_exception("wrong file format (couldn't read header)");
}

void read_header(FILE* file, int* dims)
{
	if (fscanf(file, "Matrix %d x %d\n",  &dims[0], &dims[1]) < 2)
		throw Matrix_exception("wrong file format (couldn't read header)");
}

void read_elems(int file, double* buf, int count)
{
	if (read(file, buf, count*2*sizeof(double)) < count*2*sizeof(double))
		throw Matrix_exception("wrong file format (couldn't read elements)");
}

void read_elems(FILE* file, double* buf, int count)
{
	for (int i = 0; i < count; i++)
		if (i != count-1)
		{
			if (fscanf(file, "(%lf,%lf) ", &buf[2*i], &buf[2*i+1]) < 2)
				throw Matrix_exception("wrong file format (couldn't read elements)");
		}
		else
		{
			if (fscanf(file, "(%lf,%lf)", &buf[2*i], &buf[2*i+1]) < 2)
				throw Matrix_exception("wrong file format (couldn't read elements)");
		}
}

// Write:

void write_header(ofstream& fd, int rows, int cols)
{
	fd.write("MATRIX", 6);
	fd.write(reinterpret_cast<const char*>(&rows), sizeof(int));
	fd.write(reinterpret_cast<const char*>(&cols), sizeof(int));
}

void write_header(FILE* file, int rows, int cols)
{
	fprintf(file, "Matrix %d x %d\n", rows, cols);
}

void write_elems(ofstream& fd, double* buf, int count)
{
	fd.write(reinterpret_cast<const char*>(buf), count*2*sizeof(double));
}

void write_elems(FILE* file, double* buf, int count)
{
	for (int i = 0; i < count; i++)
	{
		fprintf(file, "(%.3f,%.3f)", buf[2*i], buf[2*i+1]);
		if (i != count-1)
			fprintf(file, " ");
	}
}

// I/O with files:

int Matrix::readf(const char* filename, int row_block = R_BLOCK_SIZE, int col_block = C_BLOCK_SIZE)
// Read from file
// root returns 0 in case of error
// root returns 1 in case of binary file read
// root returns 2 in case of text file read
{

	FILE* file;
	int fd;
	int mode = 0;
	int dims[2];

	if (ProcessorGrid::is_root())
	{
		mode = bin_or_txt (filename);

		if (mode == 1) /* Binary file */
		{
			fd = open(filename, O_RDONLY);
			if (fd < 0)
			{
				throw Matrix_exception("failed to open file to import");
				//MPI_Abort(MPI_COMM_WORLD,-1);
			}
			read_header(fd,dims);
		}
		else

		if (mode == 2) /* Text file */
		{	
			file = fopen(filename ,"r");
			if (file == NULL)
			{
				throw Matrix_exception("failed to open file to import");
				//MPI_Abort(MPI_COMM_WORLD,-1);
			}
			read_header(file,dims);
		}
	}

	MPI_Bcast(&mode,1,MPI_INT,ProcessorGrid::root,MPI_COMM_WORLD);
	if (mode == 0)
		throw Matrix_exception("wrong file format (undefined header)");

	MPI_Bcast(dims,2,MPI_INT,ProcessorGrid::root,MPI_COMM_WORLD);
	init(dims[0],dims[1],row_block,col_block);
	int* local_proc_info = gather_info();

	if (ProcessorGrid::is_root())
	// Root
	{
		int root_row = 0;

		for (int ofs = 0; ofs < global_n_rows()*global_n_cols(); )
		{
			int row = ofs / global_n_cols();
			int col = ofs % global_n_cols();
			int trg_proc = get_proc_by_index(local_proc_info,row,col);
			// Get number of cols of local matrix on trg_proc:
			int length = local_proc_info[trg_proc*4+3];
			double* buf = (double*)malloc(length*2*sizeof(double));

			if (mode == 1) // Binary
				read_elems(fd,buf,length);
			if (mode == 2) // Text
			{
				read_elems(file,buf,length);				
				if (ofs % global_n_cols() == 0)
					fscanf(file, "\n");
				else
					fscanf(file, " ");
			}
			ofs += length;

			if (trg_proc == ProcessorGrid::root)
			{
				set_row(buf,root_row);
				root_row++;
			}
			else
			{
				MPI_Send(buf,2*length,MPI_DOUBLE,trg_proc,0,MPI_COMM_WORLD);
			}

			free(buf);
		}

		if (mode == 1)
			close(fd);
		if (mode == 2)
			fclose(file);
	}
	else
	// Not root
	{
		double* buf = (double*)malloc(n_cols*2*sizeof(double));
		MPI_Status status;
		for (int i = 0; i < n_rows; i++)
		{
			MPI_Recv(buf,n_cols*2,MPI_DOUBLE,ProcessorGrid::root,0,MPI_COMM_WORLD,&status);
			set_row(buf,i);
		}
		free(buf);
	}

	//local_data_transpose();

	return mode;
}

void Matrix::writef(int mode, const char* filename)
// Write to file:
// Variable 'mode':
// 1 in case of binary file write
// 2 in case of text file write
{
	if (mode != 1 && mode != 2)
		throw Matrix_exception("undefined export mode");

	ofstream fd;
	FILE* file;

	//local_data_transpose();

	if (ProcessorGrid::is_root())
	{
		if (mode == 1) /* Binary file */
		{
			fd.open(filename, ios::out | ios::binary | ios::trunc);
			if (!fd.is_open())
			{
				throw Matrix_exception("failed to open file to export");
				//MPI_Abort(MPI_COMM_WORLD,-1);
			}
			write_header(fd,global_n_rows(),global_n_cols());
		}
		else

		if (mode == 2) /* Text file */
		{	
			file = fopen(filename ,"w");
			if (file == NULL)
			{
				throw Matrix_exception("failed to open file to export");
				//MPI_Abort(MPI_COMM_WORLD,-1);
			}
			write_header(file,global_n_rows(),global_n_cols());
		}
	}

	int* local_proc_info = gather_info();

	if (ProcessorGrid::is_root())
	// Root
	{
		int root_row = 0;

		for (int ofs = 0; ofs < global_n_rows()*global_n_cols(); )
		{
			int row = ofs / global_n_cols();
			int col = ofs % global_n_cols();
			int src_proc = get_proc_by_index(local_proc_info,row,col);
			int length = local_proc_info[src_proc*4+3];
			double* buf = (double*)malloc(length*2*sizeof(double));

			if (src_proc == ProcessorGrid::root)
			{
				get_row(buf,root_row);
				root_row++;
			}
			else
			{
				MPI_Status status;
				MPI_Recv(buf,2*length,MPI_DOUBLE,src_proc,0,MPI_COMM_WORLD,&status);
			}
			ofs += length;

			if (mode == 1) // Binary
				write_elems(fd,buf,length);
			if (mode == 2) // Text
			{
				write_elems(file,buf,length);
				if (ofs % global_n_cols() == 0)
					fprintf(file, "\n");
				else
					fprintf(file, " ");
			}

			free(buf);
		}

		if (mode == 1)
			fd.close();
		if (mode == 2)
			fclose(file);
	}
	else
	// Not root
	{
		double* buf = (double*)malloc(n_cols*2*sizeof(double));
		for (int i = 0; i < n_rows; i++)
		{
			get_row(buf,i);
			MPI_Send(buf,n_cols*2,MPI_DOUBLE,ProcessorGrid::root,0,MPI_COMM_WORLD);
		}
		free(buf);
	}
}

// Stream conditional print

void Matrix::print_on_condition(ostream& out, bool (*condition)(int i, int j))
{
	int min_dim = (global_n_rows() < global_n_cols())? global_n_rows() : global_n_cols();
	out << flush;
	ProcessorGrid::barrier();
	for (int i = 0; i < global_n_rows(); i++)
		for (int j = 0; j < global_n_cols(); j++)
		{
			if (in_block(i,j) && condition(i,j))
			{
				out << data[(i - info.row_offset()) + n_rows * (j - info.col_offset())];
				if (j == (global_n_cols() - 1))
					out << endl;
				else
					out << ' ';
				out << flush;
			}
			ProcessorGrid::barrier();
		}
}

bool daigonal_elements(int i, int j)
{
	return (i == j);
}

bool all_elements(int i, int j)
{
	return true;
}

void Matrix::print_diagonal_abs(FILE* file)
// Need serialization here
{
	int min_dim = (global_n_rows() < global_n_cols())? global_n_rows() : global_n_cols();
	int* local_proc_info = gather_info();

	for (int i = 0; i < min_dim; i++)
		if (ProcessorGrid::is_root())
		{
			int src_proc = get_proc_by_index(local_proc_info,i,i);
			double magnitude;
			if (src_proc != ProcessorGrid::root)
			{
				MPI_Status status;
				MPI_Recv(&magnitude,1,MPI_DOUBLE,src_proc,0,MPI_COMM_WORLD,&status);
			}
			else
			{
				magnitude = abs(data[(i - info.row_offset()) + n_rows * (i - info.col_offset())]);
			}

			fprintf(file, "%lf", magnitude);
			if (i != min_dim-1)
				fprintf(file, " ");
			else
				fprintf(file, "\n");
		}
		else
			if (in_block(i,i))
			{
				double magnitude = abs(data[(i - info.row_offset()) + n_rows * (i - info.col_offset())]);
				MPI_Request request;
				MPI_Isend(&magnitude,1,MPI_DOUBLE,ProcessorGrid::root,0,MPI_COMM_WORLD,&request);
			}
}

// Matrix generator

void Matrix::generate (complexd (*func)(int i, int j))
{
	for (int i = 0; i < global_n_rows(); i++)
		for (int j = 0; j < global_n_cols(); j++)
			if (in_block(i,j))
				data[(i - info.row_offset()) + n_rows * (j - info.col_offset())] = func(i,j);
}

double normalized_rand()
{
	return static_cast <double> (rand()) / static_cast <double> (RAND_MAX);
}

complexd c_normalized_rand()
{
	return complexd(normalized_rand(),normalized_rand());
}

complexd identity_matrix(int i, int j)
{
	return (i == j)? 1 : 0;
}

complexd diagonal_natural_sequence(int i, int j)
{
	return (i == j)? i : 0;
}

complexd zero_matrix(int i, int j)
{
	return 0;
}

complexd random_lower_triangle(int i, int j)
{
	if (i >= j)
		return c_normalized_rand();
	else
		return 0;
}

complexd index_indicator(int i, int j)
{
	return complexd(i,j);
}

complexd test_H(int i, int j)
{
	if (i == 0 && j == 1)
		return 1;
	if (i == 1 && j == 0)
		return 1;
	if (i == j)
		return 1;
	return 0;
}

complexd test_R(int i, int j)
{
	return (i == 0 && j == 0)? 1 : 0;
}
