//view plain copy

/***********************Contents****************************************/
/*PrincipalComponentsAnalysis:C, 638lines.****************************/
/*Sampleinputdataset(final36lines).*********************************/
/**************************************************************************/

/*********************************/
/*PrincipalComponentsAnalysis*/
/*********************************/

/*********************************************************************/
/*PrincipalComponentsAnalysisortheKarhunen-Loeveexpansionisa
classicalmethodfordimensionalityreductionorexploratorydata
analysis.Onereferenceamongmanyis:F.MurtaghandA.Heck,
MultivariateDataAnalysis,KluwerAcademic,Dordrecht,1987.*/

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <stdlib.h>

#define SIGN(a, b)((b)<0?-fabs(a):fabs(a))//根据b的符号确定a的正负

main(int argc, char *argv[])
{
	FILE *stream;
	int n, m, i, j, k, k2;
	double **symmat2;
	float **data, **matrix(), **symmat, *vector(), *evals, *interm;
	void free_matrix(), free_vector(), corcol(), covcol(), scpcol();
	void tred2(), tqli();
	float in_value;
	char option, *strncpy();

	/*********************************************************************
	Getfromcommandline:
	   inputdatafilename,#rows,#cols,option.
		  
		   Openinputfile:fopenopensthefilewhosenameisstoredinthe
			  pointerargv[argc-1];ifunsuccessful,errormessageisprintedto
				 stderr.
					*********************************************************************/

	if(argc!= 5)
	{
		printf("Syntaxhelp:PCAfilename#rows#colsoption\n\n");
		printf("(filename--givefullpathname,\n");
		printf("#rows\n");
		printf("#cols--integervalues,\n");
		printf("option--R(recommended)forcorrelationanalysis,\n");
		printf("Vforvariance/covarianceanalysis\n");
		printf("SforSSCPanalysis.)\n");
		exit(1);
	}

	n= atoi(argv[2]);/*#rows*/
	m= atoi(argv[3]);/*#columns*/
	strncpy(&option, argv[4], 1);/*Analysisoption*/

	printf("No.ofrows:%d,no.ofcolumns:%d.\n", n, m);
	printf("Inputfile:%s.\n", argv[1]);

	if((stream= fopen(argv[1], "r")) == NULL)
	{
		fprintf(stderr, "Program%s:cannotopenfile%s\n",
			argv[0], argv[1]);
		fprintf(stderr, "Exitingtosystem.");
		exit(1);
		/*Note:inversionsofDOSpriorto3.0,argv[0]containsthe
			  string"C".*/
	}

	/*Nowreadindata.*/

	data= matrix(n, m);/*Storageallocationforinputdata*/

	for(i= 1;i<= n;i++)
	{
		for(j= 1;j<= m;j++)
		{
			fscanf(stream, "%f", &in_value);
			data[i][j] = in_value;
		}
	}

	/*Checkon(partof)inputdata.
		for(i=1;i<=18;i++){for(j=1;j<=8;j++){
			printf("%7.1f",data[i][j]);}printf("\n");}
					*/

	symmat= matrix(m, m);/*Allocationofcorrelation(etc.)matrix*/

							 /*Lookatanalysisoption;branchinaccordancewiththis.*/

	switch (option)
	{
		case'R':
		case'r' :
			printf("Analysisofcorrelationschosen.\n");
		corcol(data, n, m, symmat);

		/*Outputcorrelationmatrix.
								  for(i=1;i<=m;i++){
															for(j=1;j<=8;j++){
																					   printf("%7.4f",symmat[i][j]);}
																												   printf("\n");}
																																			   */
		break;
		case'V':
		case'v' :
			printf("Analysisofvariances-covarianceschosen.\n");
		covcol(data, n, m, symmat);

		/*Outputvariance-covariancematrix.
								  for(i=1;i<=m;i++){
															for(j=1;j<=8;j++){
																					  printf("%7.1f",symmat[i][j]);}
																												  printf("\n");}
																																			  */
		break;
		case'S':
		case's' :
			printf("Analysisofsums-of-squares-cross-products");
		printf("matrixchosen.\n");
		scpcol(data, n, m, symmat);

		/*OutputSSCPmatrix.
								 for(i=1;i<=m;i++){
														  for(j=1;j<=8;j++){
																					printf("%7.1f",symmat[i][j]);}
																												printf("\n");}
																																			*/
		break;
	default:
		printf("Option:%s\n", option);
		printf("Foroption,pleasetypeR,V,orS\n");
		printf("(upperorlowercase).\n");
		printf("Exitingtosystem.\n");
		exit(1);
		break;
	}

	/*********************************************************************
	Eigen-reduction
		**********************************************************************/

		/*Allocatestoragefordummyandnewvectors.*/
	evals= vector(m);/*Storagealloc.forvectorofeigenvalues*/
	interm= vector(m);/*Storagealloc.for'intermediate'vector*/
	symmat2= matrix(m, m);/*Duplicateofcorrelation(etc.)matrix*/
	for(i= 1;i<= m;i++) {
		for(j= 1;j<= m;j++) {
			symmat2[i][j] = symmat[i][j];/*Neededbelowforcol.projections*/
		}
	}
	tred2(symmat, m, evals, interm);/*Triangulardecomposition*/
	tqli(evals, interm, m, symmat);/*Reductionofsym.trid.matrix*/
									  /*evalsnowcontainstheeigenvalues,
										  columnsofsymmatnowcontaintheassociatedeigenvectors.*/

	printf("\nEigenvalues:\n");
	for(j= m;j>= 1;j--) {
		printf("%18.5f\n", evals[j]);
	}
	printf("\n(Eigenvaluesshouldbestrictlypositive;limited\n");
	printf("precisionmachinearithmeticmayaffectthis.\n");
	printf("Eigenvaluesareoftenexpressedascumulative\n");
	printf("percentages,representingthe'percentagevariance\n");
	printf("explained'bytheassociatedaxisorprincipalcomponent.)\n");

	printf("\nEigenvectors:\n");
	printf("(Firstthree;theirdefinitionintermsoforiginalvbes.)\n");
	for(j= 1;j<= m;j++) {
		for(i= 1;i<= 3;i++) {
			printf("%12.4f", symmat[j][m - i + 1]);
		}
		printf("\n");
	}

	/*Formprojectionsofrow-pointsonfirstthreeprin.components.*/
		 /*Storein'data',overwritingoriginaldata.*/
	for(i= 1;i<= n;i++) {
		for(j= 1;j<= m;j++) {
			interm[j] = data[i][j];
		}/*data[i][j]willbeoverwritten*/
		for(k= 1;k<= 3;k++) {
			data[i][k] = 0.0;
			for(k2= 1;k2<= m;k2++) {
				data[i][k] += interm[k2] * symmat[k2][m - k + 1];
			}
		}
	}

	printf("\nProjectionsofrow-pointsonfirst3prin.comps.:\n");
	for(i= 1;i<= n;i++) {
		for(j= 1;j<= 3;j++) {
			printf("%12.4f", data[i][j]);
		}
		printf("\n");
	}

	/*Formprojectionsofcol.-pointsonfirstthreeprin.components.*/
		 /*Storein'symmat2',overwritingwhatwasstoredinthis.*/
	for(j= 1;j<= m;j++) {
		for(k= 1;k<= m;k++) {
			interm[k] = symmat2[j][k];
		}/*symmat2[j][k]willbeoverwritten*/
		for(i= 1;i<= 3;i++) {
			symmat2[j][i] = 0.0;
			for(k2= 1;k2<= m;k2++) {
				symmat2[j][i] += interm[k2] * symmat[k2][m - i + 1];
			}
			if(evals[m - i + 1]>0.0005)/*Guardagainstzeroeigenvalue*/
				symmat2[j][i] /= sqrt(evals[m - i + 1]);/*Rescale*/
			else
				symmat2[j][i] = 0.0;/*Standardkludge*/
		}
	}

	printf("\nProjectionsofcolumn-pointsonfirst3prin.comps.:\n");
	for(j= 1;j<= m;j++) {
		for(k= 1;k<= 3;k++) {
			printf("%12.4f", symmat2[j][k]);
		}
		printf("\n");
	}

	free_matrix(data, n, m);
	free_matrix(symmat, m, m);
	free_matrix(symmat2, m, m);
	free_vector(evals, m);
	free_vector(interm, m);

}

/**Correlationmatrix:creation***********************************/

void corcol(data, n, m, symmat)
float**data, **symmat;
int n, m;
/*Createm*mcorrelationmatrixfromgivenn*mdatamatrix.*/
{
	float eps= 0.005;
	double x, *stddev;
	float  *mean,  *vector();
	int i, j, j1, j2;

	/*Allocatestorageformeanandstd.dev.vectors*/

	mean= vector(m);
	stddev= vector(m);

	/*Determinemeanofcolumnvectorsofinputdatamatrix*/

	for(j= 1;j<= m;j++)
	{
		mean[j] = 0.0;
		for(i= 1;i<= n;i++)
		{
			mean[j] += data[i][j];
		}
		mean[j] /= (float)n;
	}

	printf("\nMeansofcolumnvectors:\n");
	for(j= 1;j<= m;j++) {
		printf("%7.1f", mean[j]);
	}printf("\n");

	/*Determinestandarddeviationsofcolumnvectorsofdatamatrix.*/

	for(j= 1;j<= m;j++)
	{
		stddev[j] = 0.0;
		for(i= 1;i<= n;i++)
		{
			stddev[j] += ((data[i][j] - mean[j])*
				(data[i][j] - mean[j]));
		}
		stddev[j] /= (float)n;
		stddev[j] = sqrt(stddev[j]);
		/*Thefollowinginaninelegantbutusualwaytohandle
				near-zerostd.dev.values,whichbelowwouldcauseazero-
						divide.*/
		if(stddev[j] <= eps)stddev[j] = 1.0;
	}

	printf("\nStandarddeviationsofcolumns:\n");
	for(j= 1;j<= m;j++) { printf("%7.1f", stddev[j]); }
	printf("\n");

	/*Centerandreducethecolumnvectors.*/

	for(i= 1;i<= n;i++)
	{
		for(j= 1;j<= m;j++)
		{
			data[i][j] -= mean[j];
			x= sqrt((float)n);
			x*= stddev[j];
			data[i][j] /= x;
		}
	}

	/*Calculatethem*mcorrelationmatrix.*/
	for(j1= 1;j1<= m - 1;j1++)
	{
		symmat[j1][j1] = 1.0;
		for(j2= j1 + 1;j2<= m;j2++)
		{
			symmat[j1][j2] = 0.0;
			for(i= 1;i<= n;i++)
			{
				symmat[j1][j2] += (data[i][j1] * data[i][j2]);
			}
			symmat[j2][j1] = symmat[j1][j2];
		}
	}
	symmat[m][m] = 1.0;

	return;

}

/**Variance-covariancematrix:creation*****************************/

void covcol(data, n, m, symmat)
float **data, **symmat;
int n, m;
/*Createm*mcovariancematrixfromgivenn*mdatamatrix.*/
{
	float *mean, *vector();
	int i, j, j1, j2;

	/*Allocatestorageformeanvector*/

	mean= vector(m);

	/*Determinemeanofcolumnvectorsofinputdatamatrix*/

	for(j= 1;j<= m;j++)
	{
		mean[j] = 0.0;
		for(i= 1;i<= n;i++)
		{
			mean[j] += data[i][j];
		}
		mean[j] /= (float)n;
	}

	printf("\nMeansofcolumnvectors:\n");
	for(j= 1;j<= m;j++) {
		printf("%7.1f", mean[j]);
	}printf("\n");

	/*Centerthecolumnvectors.*/

	for(i= 1;i<= n;i++)
	{
		for(j= 1;j<= m;j++)
		{
			data[i][j] -= mean[j];
		}
	}

	/*Calculatethem*mcovariancematrix.*/
	for(j1= 1;j1<= m;j1++)
	{
		for(j2= j1;j2<= m;j2++)
		{
			symmat[j1][j2] = 0.0;
			for(i= 1;i<= n;i++)
			{
				symmat[j1][j2] += data[i][j1] * data[i][j2];
			}
			symmat[j2][j1] = symmat[j1][j2];
		}
	}

	return;

}

/**Sums-of-squares-and-cross-productsmatrix:creation**************/

void scpcol(data, n, m, symmat)
float**data, **symmat;
int n, m;
/*Createm*msums-of-cross-productsmatrixfromn*mdatamatrix.*/
{
	int i, j1, j2;

	/*Calculatethem*msums-of-squares-and-cross-productsmatrix.*/

	for(j1= 1;j1<= m;j1++)
	{
		for(j2= j1;j2<= m;j2++)
		{
			symmat[j1][j2] = 0.0;
			for(i= 1;i<= n;i++)
			{
				symmat[j1][j2] += data[i][j1] * data[i][j2];
			}
			symmat[j2][j1] = symmat[j1][j2];
		}
	}

	return;

}

/**Errorhandler**************************************************/

void erhand(err_msg)
char err_msg[];
/*Errorhandler*/
{
	fprintf(stderr, "Run-timeerror:\n");
	fprintf(stderr, "%s\n", err_msg);
	fprintf(stderr, "Exitingtosystem.\n");
	exit(1);
}

/**Allocationofvectorstorage***********************************/

float *vector(n)
int n;
/*Allocatesafloatvectorwithrange[1..n].*/
{

	float*v;

	v= (float*)malloc((unsigned)n * sizeof(float));
	if(!v)erhand("Allocationfailureinvector().");
	return v - 1;

}

/**Allocationoffloatmatrixstorage*****************************/

float **matrix(n, m)
int n, m;
/*Allocateafloatmatrixwithrange[1..n][1..m].*/
{
	int i;
	float **mat;

	/*Allocatepointerstorows.*/
	mat= (float**)malloc((unsigned)(n) * sizeof(float*));
	if(!mat)erhand("Allocationfailure1inmatrix().");
	mat-= 1;

	/*Allocaterowsandsetpointerstothem.*/
	for(i= 1;i<= n;i++)
	{
		mat[i] = (float*)malloc((unsigned)(m) * sizeof(float));
		if(!mat[i])erhand("Allocationfailure2inmatrix().");
		mat[i] -= 1;
	}

	/*Returnpointertoarrayofpointerstorows.*/
	return mat;

}

/**Deallocatevectorstorage*********************************/

void free_vector(v, n)
float*v;
int n;
/*Freeafloatvectorallocatedbyvector().*/
{
	free((char*)(v + 1));
}

/**Deallocatefloatmatrixstorage***************************/

void free_matrix(mat, n, m)
float **mat;
int n, m;
/*Freeafloatmatrixallocatedbymatrix().*/
{
	int i;

	for(i= n;i>= 1;i--)
	{
		free((char*)(mat[i] + 1));
	}
	free((char*)(mat + 1));
}

/**Reduceareal,symmetricmatrixtoasymmetric,tridiag.matrix.*/

void tred2(a, n, d, e)
float**a, *d, *e;
/*float**a,d[],e[];*/
int n;
/*Householderreductionofmatrixatotridiagonalform.
Algorithm:Martinetal.,Num.Math.11,181-195,1968.
   Ref:Smithetal.,MatrixEigensystemRoutines--EISPACKGuide
	  Springer-Verlag,1976,pp.489-494.
			  WHPressetal.,NumericalRecipesinC,CambridgeUP,
					  1988,pp.373-374.*/
{
	int l, k, j, i;
	double g, scale;
	float  hh, h,  f;

	for(i= n;i>= 2;i--)
	{
		l= i- 1;
		h= scale= 0.0;
		if(l>1)
		{
			for(k= 1;k<= l;k++)
				scale+= fabs(a[i][k]);
			if(scale== 0.0)
				e[i] = a[i][l];
			else
			{
				for(k= 1;k<= l;k++)
			{
				a[i][k] /= scale;
				h+= a[i][k] * a[i][k];
			}
			f= a[i][l];
 			g= f>0?-sqrt(h):sqrt(h);
			e[i] = scale*g;
			h-= f*g;
			a[i][l] = f- g;
			f= 0.0;
			for(j= 1;j<= l;j++)
			{
				a[j][i] = a[i][j] / h;
				g= 0.0;
				for(k= 1;k<= j;k++)
					g+= a[j][k] * a[i][k];
				for(k= j + 1;k<= l;k++)
					g+= a[k][j] * a[i][k];
				e[j] = g/ h;
				f+= e[j] * a[i][j];
			}
			hh= f/ (h+ h);
			for(j= 1;j<= l;j++)
			{
				f= a[i][j];
				e[j] = g= e[j] - hh*f;
				for(k= 1;k<= j;k++)
					a[j][k] -= (f*e[k] + g*a[i][k]);
			}
			}
		}
		else
			e[i] = a[i][l];
		d[i] = h;
	}
	d[1] = 0.0;
	e[1] = 0.0;
	for(i= 1;i<= n;i++)
	{
		l= i- 1;
		if(d[i])
		{
			for(j= 1;j<= l;j++)
			{
				g= 0.0;
				for(k= 1;k<= l;k++)
					g+= a[i][k] * a[k][j];
				for(k= 1;k<= l;k++)
					a[k][j] -= g*a[k][i];
			}
		}
		d[i] = a[i][i];
		a[i][i] = 1.0;
		for(j= 1;j<= l;j++)
			a[j][i] = a[i][j] = 0.0;
	}
}

/**TridiagonalQLalgorithm--Implicit**********************/

void tqli(d, e, n, z)
float d[], e[], **z;
int n;
{
	int m, l, iter, i, k;
	double r, g, dd;
	float s, p,  f,  c, b;
	void erhand();

	for(i= 2;i<= n;i++)
		e[i - 1] = e[i];
	e[n] = 0.0;
	for(l= 1;l<= n;l++)
	{
		iter= 0;
		do
		{
			for(m= l;m<= n - 1;m++)
		{
			dd= fabs(d[m]) + fabs(d[m + 1]);
			if(fabs(e[m]) + dd== dd)break;
		}
		if(m!= l)
		{
			if(iter++ == 30)erhand("NoconvergenceinTLQI.");
			g= (d[l + 1] - d[l]) / (2.0*e[l]);
			r= sqrt((g*g) + 1.0);
			g= d[m] - d[l] + e[l] / (g+ SIGN(r,g));
			s= c= 1.0;
			p= 0.0;
			for(i= m - 1;i>= l;i--)
			{
				f= s*e[i];
				b= c*e[i];
				if(fabs(f) >= fabs(g))
				{
					c= g/ f;
					r= sqrt((c*c) + 1.0);
					e[i + 1] = f*r;
					c*= (s= 1.0 / r);
				}
				else
				{
					s= f/ g;
				r= sqrt((s*s) + 1.0);
				e[i + 1] = g*r;
				s*= (c= 1.0 / r);
				}
				g= d[i + 1] - p;
				r= (d[i] - g)*s+ 2.0*c*b;
				p= s*r;
				d[i + 1] = g+ p;
				g= c*r- b;
				for(k= 1;k<= n;k++)
				{
					f= z[k][i + 1];
					z[k][i + 1] = s*z[k][i] + c*f;
					z[k][i] = c*z[k][i] - s*f;
				}
			}
			d[l] = d[l] - p;
			e[l] = g;
			e[m] = 0.0;
		}
		}while(m!= l);
	}
}