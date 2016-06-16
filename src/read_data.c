/* read_data.h */
/*********************************************/
/*	read data ( dim of matrix )					*/
/*	read data ( data )								*/
/*********************************************/

#include <stdio.h>
#include <stdlib.h>
#include "read_data.h"

void ReadFail(	int status
					)
{
	if( !status ){
		printf("Reading failed.\n");
		exit(0);
	}
}

void ReadDim(	const char	*filename,
					int			*n,			/* number of data */
					int			*p,			/* number of explanatory variable */ 
					int			*i_ex			/* index of explained varaible */
					)
{
	FILE	*file;
	int	status;
	char	s[10];
	
	/* open file */
	file = fopen(filename, "r");
	if( file==NULL ){
		printf("Could not open file <%s>.\n", filename);
		exit(0);
	}

	/* skip one line */
	if( fgets( s, 10, file)==NULL ){
		printf("Error reading file <%s>.\n", filename);
	  	exit(0);
   }

	/* read n */ 
	status = fscanf( file, "%d", n);
	ReadFail(status);

	/* read p */
	status = fscanf( file, "%d", p);
	ReadFail(status);
	
	/* read i_ex */
	status = fscanf( file, "%d", i_ex);
	ReadFail(status);

	fclose( file);
}

void ReadData(	const char	*filename,
					int			n,
					int			p,
					double		*data
					)
{
	int 	i;

	FILE	*file;
	int	buf;
	int	status;
	char	s[10];
	
	/* open file */
	file = fopen(filename, "r");
	if( file==NULL ){
		printf("Could not open file <%s>.\n", filename);
		exit(0);
	}

	/* skip one line */
	if( fgets( s, 10, file)==NULL ){
		printf("Error reading file <%s>.\n", filename);
	  	exit(0);
   }

	/* skip n */
	status = fscanf( file, "%d", &buf);
	ReadFail(status);
	/* skip p */
	status = fscanf( file, "%d", &buf);
	ReadFail(status);
	/* skip i_ex */
	status = fscanf( file, "%d", &buf);
	ReadFail(status);

	/* read data */
	for(i=0; i<(n*(p+1)); ++i){
		status = fscanf( file, "%lf", (data+i));
		ReadFail(status);
		
	}

	fclose( file);

}
