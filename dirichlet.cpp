/*  Blocks are decomposed like that  *\
      proc0 \/
    ******************
    *       *        *
    * proc0 * proc1  *
    *       *        *
pr2>****************** <-proc1
    *       *        *
    * proc2 * proc3  *
    *       *        *
    ******************
            ^
            | proc3
\*                                   */


#include<stdlib.h>
#include<iostream>
#include<fstream>
#include<math.h>
#include<mpi.h>
using namespace std;
double **f_old,**f_new, step,*f_bound,center,res,error;
int nNodes,nProcs,id;
double F(double x, double y)
{
    return -2.*sin(x+y);
}
double F_bound(double x,double y)
{
    return sin(x+y);
}
double function(int i, int j)
{
    double x,y;
    if(id==0)
    {
        x=step*(i);
        y=step*(j+nNodes-1);
    }
    if(id==1)
    {
        x=step*(i+nNodes-1);
        y=step*(j+nNodes-1);
    }
    if(id==2)
    {
        x=step*(i);
        y=step*(j);
    }
    if(id==3)
    {
        x=step*(i+nNodes-1);
        y=step*(j);
    }
    return F(x,y);
}
double BC(int i, int j)
{
    double x,y;
    if(id==0)
    {
        x=step*(i);
        y=step*(j+nNodes-1);
    }
    if(id==1)
    {
        x=step*(i+nNodes-1);
        y=step*(j+nNodes-1);
    }
    if(id==2)
    {
        x=step*(i);
        y=step*(j);
    }
    if(id==3)
    {
        x=step*(i+nNodes-1);
        y=step*(j);
    }
    return F_bound(x,y);
}
void Initialize(double Square_size)
{
    int i,j;
    step=Square_size/(nNodes*2-2.);
    f_old=(double**)malloc(nNodes*sizeof(double*));
    f_new=(double**)malloc(nNodes*sizeof(double*));
    f_bound=(double*)malloc(nNodes*sizeof(double));
    for(i=0;i<nNodes;i++)
    {
        f_old[i]=(double*)malloc(nNodes*sizeof(double));
        f_new[i]=(double*)malloc(nNodes*sizeof(double));
    }
    for(i=0;i<nNodes;i++)
    {
        for(j=0;j<nNodes;j++)
        {
            f_old[i][j]=BC(i,j);
            f_new[i][j]=BC(i,j);
        }
        f_bound[i]=0.;
    }
}
double Current(int i, int j)
{
    return 1./4.*(f_old[i-1][j]+f_old[i+1][j]+f_old[i][j-1]+f_old[i][j+1]-step*step*function(i,j));
}



void boundaries_init()
{
    int i,j;
    if(id==0)
    {
        for(i=0;i<nNodes;i++)
        {
            f_old[0][i]=BC(0,i);
            f_old[i][nNodes-1]=BC(i,nNodes-1);
        }
    }
    if(id==1)
    {
        for(i=0;i<nNodes;i++)
        {
            f_old[i][nNodes-1]=BC(i,nNodes-1);
            f_old[nNodes-1][i]=BC(nNodes-1,i);
        }
    }
    if(id==2)
    {
        for(i=0;i<nNodes;i++)
        {
            f_old[i][0]=BC(i,0);
            f_old[0][i]=BC(0,i);
        }
    }
    if(id==3)
    {
        for(i=0;i<nNodes;i++)
        {
            f_old[nNodes-1][i]=BC(nNodes-1,i);
            f_old[i][0]=BC(i,0);
        }
    }
}
void boundaries()
{
    int i,j;
    if(id==0)
    {
        for(i=0;i<nNodes;i++)
        {
            f_new[0][i]=BC(0,i);
            f_new[i][nNodes-1]=BC(i,nNodes-1);
        }
    }
    if(id==1)
    {
        for(i=0;i<nNodes;i++)
        {
            f_new[i][nNodes-1]=BC(i,nNodes-1);
            f_new[nNodes-1][i]=BC(nNodes-1,i);
        }
    }
    if(id==2)
    {
        for(i=0;i<nNodes;i++)
        {
            f_new[i][0]=BC(i,0);
            f_new[0][i]=BC(0,i);
        }
    }
    if(id==3)
    {
        for(i=0;i<nNodes;i++)
        {
            f_new[nNodes-1][i]=BC(nNodes-1,i);
            f_new[i][0]=BC(i,0);
        }
    }
}

void proc_bound()
{
    int i,j;
    MPI_Status status;
    if(id==0)
    {
        for(i=0;i<nNodes;i++)
            f_bound[i]=f_old[i][1];
        MPI_Send(f_bound,nNodes,MPI_DOUBLE,2,0,MPI_COMM_WORLD);
    }
    if(id==2)
    {
        MPI_Recv(f_bound,nNodes,MPI_DOUBLE,0,0,MPI_COMM_WORLD,&status);
        for(i=1;i<nNodes-1;i++)
            f_new[i][nNodes-1]=1./4.*(f_old[i-1][nNodes-1]+f_old[i+1][nNodes-1]+f_old[i][nNodes-2]+f_bound[i]-step*step*function(i,nNodes-1));
    }
    if(id==1)
    {
        for(i=0;i<nNodes;i++)
            f_bound[i]=f_old[1][i];
        MPI_Send(f_bound,nNodes,MPI_DOUBLE,0,1,MPI_COMM_WORLD);
    }
    if(id==0)
    {
        MPI_Recv(f_bound,nNodes,MPI_DOUBLE,1,1,MPI_COMM_WORLD,&status);
        for(i=1;i<nNodes-1;i++)
            f_new[nNodes-1][i]=1./4.*(f_old[nNodes-2][i]+f_bound[i]+f_old[nNodes-1][i-1]+f_old[nNodes-1][i+1]-step*step*function(nNodes-1,i));
    }
    if(id==2)
    {
        for(i=0;i<nNodes;i++)
            f_bound[i]=f_old[nNodes-2][i];
        MPI_Send(f_bound,nNodes,MPI_DOUBLE,3,2,MPI_COMM_WORLD);
    }
    if(id==3)
    {
        MPI_Recv(f_bound,nNodes,MPI_DOUBLE,2,2,MPI_COMM_WORLD,&status);
        for(i=1;i<nNodes-1;i++)
            f_new[0][i]=1./4.*(f_bound[i]+f_old[1][i]+f_old[0][i-1]+f_old[0][i+1]-step*step*function(0,i));
    }
    if(id==3)
    {
        for(i=0;i<nNodes;i++)
            f_bound[i]=f_old[i][nNodes-2];
        MPI_Send(f_bound,nNodes,MPI_DOUBLE,1,3,MPI_COMM_WORLD);
    }
    if(id==1)
    {
        MPI_Recv(f_bound,nNodes,MPI_DOUBLE,3,3,MPI_COMM_WORLD,&status);
        for(i=1;i<nNodes-1;i++)
            f_new[i][0]=1./4.*(f_old[i-1][0]+f_old[i+1][0]+f_bound[i]+f_old[i][1]-step*step*function(i,0));
    }
    double t0=0.,t1=0.,t2=0.,t3=0.,r1=0.,r2=0.,r3=0.;
    if(id==1)
        t1=f_old[1][0];
    else
        t1=0.;
    if(id==2)
        t2=f_old[nNodes-2][nNodes-1];
    else
        t2=0.;
    if(id==3)
        t3=f_old[0][nNodes-2];
    else
        t3=0.;
    MPI_Reduce(&t1,&r1,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
    MPI_Reduce(&t2,&r2,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
    MPI_Reduce(&t3,&r3,1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
    if(id==0)
    {
        t0=f_old[nNodes-1][1];
        center=1./4.*(t0+r1+r2+r3-step*step*function(nNodes-1,0));
        f_new[nNodes-1][0]=center;
    }
    MPI_Bcast(&center,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
    if(id==1)
        f_new[0][0]=center;
    if(id==2)
        f_new[nNodes-1][nNodes-1]=center;
    if(id==3)
        f_new[0][nNodes-1]=center;
}

void bcast()
{
    int i,j;
    MPI_Status status;
    if(id==0)
    {
        for(i=0;i<nNodes;i++)
            f_bound[i]=f_new[nNodes-1][i];
        MPI_Send(f_bound,nNodes,MPI_DOUBLE,1,5,MPI_COMM_WORLD);
    }
    if(id==1)
    {
        MPI_Recv(f_bound,nNodes,MPI_DOUBLE,0,5,MPI_COMM_WORLD,&status);
        for(i=1;i<nNodes-1;i++)
            f_new[0][i]=f_bound[i];
    }
    if(id==1)
    {
        for(i=0;i<nNodes;i++)
            f_bound[i]=f_new[i][0];
        MPI_Send(f_bound,nNodes,MPI_DOUBLE,3,6,MPI_COMM_WORLD);
    }
    if(id==3)
    {
        MPI_Recv(f_bound,nNodes,MPI_DOUBLE,1,6,MPI_COMM_WORLD,&status);
        for(i=1;i<nNodes-1;i++)
            f_new[i][nNodes-1]=f_bound[i];
    }
    if(id==3)
    {
        for(i=0;i<nNodes;i++)
            f_bound[i]=f_new[0][i];
        MPI_Send(f_bound,nNodes,MPI_DOUBLE,2,7,MPI_COMM_WORLD);
    }
    if(id==2)
    {
        MPI_Recv(f_bound,nNodes,MPI_DOUBLE,3,7,MPI_COMM_WORLD,&status);
        for(i=1;i<nNodes-1;i++)
            f_new[nNodes-1][i]=f_bound[i];
    }
    if(id==2)
    {
        for(i=0;i<nNodes;i++)
            f_bound[i]=f_new[i][nNodes-1];
        MPI_Send(f_bound,nNodes,MPI_DOUBLE,0,8,MPI_COMM_WORLD);
    }
    if(id==0)
    {
        MPI_Recv(f_bound,nNodes,MPI_DOUBLE,2,8,MPI_COMM_WORLD,&status);
        for(i=1;i<nNodes-1;i++)
            f_new[i][0]=f_bound[i];
    }
}

void residuals()
{
    int i,j;
    double Residual=0.,Error=0.,x,y;
    for(i=1;i<nNodes-1;i++)
    {
        for(j=1;j<nNodes-1;j++)
        {
            if(Residual<fabs(f_new[i][j]-f_old[i][j]))
                Residual=fabs(f_new[i][j]-f_old[i][j]);
            if(Error<fabs(f_new[i][j]-BC(i,j)))
                Error=fabs(f_new[i][j]-BC(i,j));
        }
    }
    res=0.;
    error=0.;
    MPI_Allreduce(&Residual,&res,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    MPI_Allreduce(&Error,&error,1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
    res=res/4.;
    error=error/4.;
    if(id==0)
    {
        for(i=1;i<nNodes;i++)
        {
            for(j=0;j<nNodes-1;j++)
                f_old[i][j]=f_new[i][j];
        }
    }
    if(id==1)
    {
        for(i=0;i<nNodes-1;i++)
        {
            for(j=0;j<nNodes-1;j++)
                f_old[i][j]=f_new[i][j];
        }
    }
    if(id==2)
    {
        for(i=1;i<nNodes;i++)
        {
            for(j=1;j<nNodes;j++)
                f_old[i][j]=f_new[i][j];
        }
    }
    if(id==3)
    {
        for(i=0;i<nNodes-1;i++)
        {
            for(j=1;j<nNodes;j++)
                f_old[i][j]=f_new[i][j];
        }
    }
}

void write()
{
    int i,j;
    double A[nNodes*2-1][nNodes*2-1],Atemp[nNodes*2-1][nNodes*2-1];
    for(i=0;i<nNodes*2-1;i++)
    {
        for(j=0;j<nNodes*2-1;j++)
        {
            A[i][j]=0.;
            Atemp[i][j]=0.;
        }
    }
    for(i=0;i<nNodes;i++)
    {
        for(j=0;j<nNodes;j++)
        {
            if(id==0)
            {
                Atemp[i][j+nNodes-1]=f_new[i][j];
            }
            if((id==1)&&(i!=0))
            {
                Atemp[nNodes-1+i][j+nNodes-1]=f_new[i][j];
            }
            if((id==2)&&(i!=nNodes-1)&&(j!=nNodes-1))
            {
                Atemp[i][j]=f_new[i][j];
            }
            if((id==3)&&(j!=nNodes-1))
            {
                Atemp[i+nNodes-1][j]=f_new[i][j];
            }
        }
    }
    for(i=0;i<nNodes*2-1;i++)
    {
        for(j=0;j<nNodes*2-1;j++)
            MPI_Reduce(&Atemp[i][j],&A[i][j],1,MPI_DOUBLE,MPI_SUM,0,MPI_COMM_WORLD);
    }
    ofstream Tecplot;
    if(id==0)
    {
        Tecplot.open("out.plt");
        Tecplot<<"Variables=\"x\",\"y\",\"F(x,y)\",\"F<sub>Analythic</sub>(x,y)\",\"step\"\nZone i=	"<<nNodes*2-1<<"	j=	"<<nNodes*2-1;
        for(i=0;i<nNodes*2-1;i++)
        {
            for(j=0;j<nNodes*2-1;j++)
                Tecplot<<endl<<i*step<<"	"<<j*step<<"	"<<A[i][j]<<"	"<<sin(i*step+j*step)<<"	"<<step;
        }
        Tecplot.close();
    }
}
int main(int argc, char *argv[])
{
    int i,j,count;
    double Square_size = 6.28;
//    Square_size=Square_size;
    nNodes=100;                                 //Число точек
    Initialize(Square_size);
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD,&nProcs);
    MPI_Comm_rank(MPI_COMM_WORLD,&id);
    if(nProcs!=4)
    {
        if(id==0)
            cout<<"Wrong number of processors - must be 4"<<endl;
        MPI_Finalize();
        exit(0);
    }
    count=0;
    res=1.;
    boundaries_init();
    while((res>1e-15)&&(count<100000))
    {
        count++;
        if((id==0)&&(count>50)&&(count%50==0))
            cout<<"---------------------------------------------------------\n| It	|	Residuals	|	Error		|\n---------------------------------------------------------"<<endl;
        boundaries();
        for(i=1;i<nNodes-1;i++)
        {
            for(j=1;j<nNodes-1;j++)
                f_new[i][j]=Current(i,j);
        }
        proc_bound();
        bcast();
        residuals();
        if(id==0)
            cout<<"|"<<count<<"	|	"<<res<<"	|	"<<error<<"	|"<<endl;
    }
    write();
    if(id==0)
    cout<<"========================================================="<<endl;
    MPI_Finalize();
    return 0;
}