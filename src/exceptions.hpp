#ifndef __EXCEPTIONS_H__
#define __EXCEPTIONS_H__

#include <exception>
#include <cstdio>

using namespace std;

// Matrix exception class-----------------------------------------------------------------

class Matrix_exception: public std::exception
{
	mutable char* errstr; 

	public:

	Matrix_exception(const char* str = "")
	{
		errstr = const_cast <char*> (str);
	}
	~Matrix_exception() throw()
	{
		delete [] errstr;
	}
	virtual const char* what() const throw()
	{
		char* tmp = errstr;
		char* prefix = const_cast <char*> ("Complex matrix error: ");
		try
		{
			errstr = new char [strlen(prefix)+strlen(errstr)+2];
		}
		catch (std::exception& smth)
		{
			return "Couldn't generate an error message (there is no memory)\n";
		}
		sprintf(errstr, "%s%s.\n", prefix, tmp);
		return errstr;
	}	
};

// Solver exception class-----------------------------------------------------------------

class Solver_exception: public std::exception
{
	mutable char* errstr; 

	public:

	Solver_exception(const char* str = "")
	{
		errstr = const_cast <char*> (str);
	}
	~Solver_exception() throw()
	{
		delete [] errstr;
	}
	virtual const char* what() const throw()
	{
		char* tmp = errstr;
		char* prefix = const_cast <char*> ("Solver class error: ");
		try
		{
			errstr = new char [strlen(prefix)+strlen(errstr)+2];
		}
		catch (std::exception& smth)
		{
			return "Couldn't generate an error message (there is no memory)\n";
		}
		sprintf(errstr, "%s%s.\n", prefix, tmp);
		return errstr;
	}	
};

// Processors exception class-----------------------------------------------------------------

class Processors_exception: public std::exception
{
	mutable char* errstr; 

	public:

	Processors_exception(const char* str = "")
	{
		errstr = const_cast <char*> (str);
	}
	~Processors_exception() throw()
	{
		delete [] errstr;
	}
	virtual const char* what() const throw()
	{
		char* tmp = errstr;
		char* prefix = const_cast <char*> ("Error in processor grid: ");
		try
		{
			errstr = new char [strlen(prefix)+strlen(errstr)+2];
		}
		catch (std::exception& smth)
		{
			return "Couldn't generate an error message (there is no memory)\n";
		}
		sprintf(errstr, "%s%s.\n", prefix, tmp);
		return errstr;
	}	
};

#endif

// __EXCEPTIONS_H__