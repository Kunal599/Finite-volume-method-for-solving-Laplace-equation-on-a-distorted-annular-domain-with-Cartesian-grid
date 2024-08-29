
n_rad=51;
n_tan=51;
Imax=n_rad ;
Jmax=n_tan ;
N=5000; %number of iterations
RX=1;
RY=1;
RW=0.1;
RZ=0.1;
ThetA_WX=0;
ThetA_ZY=90;
RMS_acceptance=0.00001;
relaxation=1.50;

%%       conversion of polar to cartesian, Exact solution of PHI
R(1:Imax,1:Jmax)=0;
ThetA(1:Imax,1:Jmax)=0;
X(1:Imax,1:Jmax)=0;
Y(1:Imax,1:Jmax)=0;

for j=1:1:Jmax
    r=RW + (RZ-RW)*(j-1)/(Jmax-1);
    D_r= ((RX-RW)/(Imax-1))+(((RY-RZ-RX+RW)/(Imax-1))*(j-1)/(Jmax-1));
    D_theta= (ThetA_ZY-ThetA_WX)/(Jmax-1);
    for i=1:1:Imax
        R(i,j)=r+D_r*(i-1);
        ThetA(i,j)= ThetA_WX + D_theta*(j-1);

        X(i,j)=R(i,j)*cos(ThetA(i,j)*pi/180);
        Y(i,j)=R(i,j)*sin(ThetA(i,j)*pi/180);
        Exact_Soln(i,j)= sin(ThetA(i,j)*pi/180)/R(i,j); 
        PHI_Soln(i,j)=Exact_Soln(i,j);
    end
end

%%           boundary conditions of PHI

for i=1:Imax
    PHI_Soln(i,1)=0;
    PHI_Soln(i,Jmax)=1/R(i,Jmax);
end
for j=1:Jmax
    PHI_Soln(1,j)=sin(ThetA(1,j)*pi/180)/R(1,j);
    PHI_Soln(Imax,j)=sin(ThetA(Imax,j)*pi/180)/R(Imax,j);
end

%%          Solving PDE equation
M_AB(2:Imax-1,2:Jmax-1)=0;
N_AB(2:Imax-1,2:Jmax-1)=0;
M_BC(2:Imax-1,2:Jmax-1)=0;
N_BC(2:Imax-1,2:Jmax-1)=0;
M_CD(2:Imax-1,2:Jmax-1)=0;
N_CD(2:Imax-1,2:Jmax-1)=0;
M_DA(2:Imax-1,2:Jmax-1)=0;
N_DA(2:Imax-1,2:Jmax-1)=0;

for j=2:1:Jmax-1
    for i=2:1:Imax-1
        XA=0.25*(X(i,j)+X(i-1,j)+X(i-1,j-1)+X(i,j-1));
        YA=0.25*(Y(i,j)+Y(i-1,j)+Y(i-1,j-1)+Y(i,j-1));
        XB=0.25*(X(i+1,j)+X(i,j)+X(i+1,j-1)+X(i,j-1));
        YB=0.25*(Y(i+1,j)+Y(i,j)+Y(i+1,j-1)+Y(i,j-1));
        XC=0.25*(X(i+1,j+1)+X(i,j+1)+X(i+1,j)+X(i,j));
        YC=0.25*(Y(i+1,j+1)+Y(i,j+1)+Y(i+1,j)+Y(i,j));
        XD=0.25*(X(i,j+1)+X(i-1,j+1)+X(i,j)+X(i-1,j));
        YD=0.25*(Y(i,j+1)+Y(i-1,j+1)+Y(i,j)+Y(i-1,j));

        Side_AB= abs((XB-XA)*(Y(i,j)-Y(i,j-1))-(YB-YA)*(X(i,j)-X(i,j-1)));
        M_AB(i,j)=((XB-XA)*(XB-XA)+(YB-YA)*(YB-YA))/Side_AB;
        N_AB(i,j)=((XB-XA)*(X(i,j)-X(i,j-1)) + (YB-YA)*(Y(i,j)-Y(i,j-1)))/Side_AB;

        Side_BC= abs((XC-XB)*(Y(i,j)-Y(i+1,j))-(YC-YB)*(X(i,j)-X(i+1,j)));
        M_BC(i,j)=((XC-XB)*(XC-XB)+(YC-YB)*(YC-YB))/Side_BC;
        N_BC(i,j)=((XC-XB)*(X(i,j)-X(i+1,j)) + (YC-YB)*(Y(i,j)-Y(i+1,j)))/Side_BC;

        Side_CD= abs((XD-XC)*(Y(i,j)-Y(i,j+1))-(YD-YC)*(X(i,j)-X(i,j+1)));
        M_CD(i,j)=((XD-XC)*(XD-XC)+(YD-YC)*(YD-YC))/Side_CD;
        N_CD(i,j)=((XD-XC)*(X(i,j)-X(i,j+1)) + (YD-YC)*(Y(i,j)-Y(i,j+1)))/Side_CD;

        Side_DA= abs((XA-XD)*(Y(i,j)-Y(i-1,j))-(YA-YD)*(X(i,j)-X(i-1,j)));
        M_DA(i,j)=((XA-XD)*(XA-XD)+(YA-YD)*(YA-YD))/Side_DA;
        N_DA(i,j)=((XA-XD)*(X(i,j)-X(i-1,j)) + (YA-YD)*(Y(i,j)-Y(i-1,j)))/Side_DA;
    end
end

n=1;
while n<N+1
    SUM=0;    
    for j=2:1:Jmax-1
        for i=2:1:Imax-1
            PHD=(M_AB(i,j)+M_BC(i,j)+M_CD(i,j)+M_DA(i,j));
            PHD_1=0.25*(N_CD(i,j)- N_DA(i,j));
            PHD_2=M_CD(i,j)+0.25*(N_BC(i,j)- N_DA(i,j));
            PHD_3=0.25*(N_BC(i,j)- N_CD(i,j));
            PHD_4=M_BC(i,j)+0.25*(N_AB(i,j)- N_CD(i,j));
            PHD_5=0.25*(N_AB(i,j)- N_BC(i,j));
            PHD_6=M_AB(i,j)+0.25*(N_DA(i,j)- N_BC(i,j));
            PHD_7=0.25*(N_DA(i,j)- N_AB(i,j));
            PHD_8=M_DA(i,j)+0.25*(N_CD(i,j)- N_AB(i,j));
           
            DIFF=(PHD_1*PHI_Soln(i-1,j+1)+PHD_2*PHI_Soln(i,j+1)+PHD_3*PHI_Soln(i+1,j+1)+PHD_4*PHI_Soln(i+1,j)+PHD_5*PHI_Soln(i+1,j-1)+PHD_6*PHI_Soln(i,j-1)+PHD_7*PHI_Soln(i-1,j-1)+PHD_8*PHI_Soln(i-1,j))/PHD;
            
            SUM=SUM+(DIFF-PHI_Soln(i,j))*(DIFF-PHI_Soln(i,j));
            PHI_Soln(i,j)=PHI_Soln(i,j)+ relaxation*(DIFF-PHI_Soln(i,j));
        end
    end

            RMS=sqrt(SUM/((Imax-2)*(Jmax-2)));

            Iteration=n;

            if RMS<RMS_acceptance
                n=N+1;
            elseif n>N-1
                fprintf('convergence not achieved in %d steps for %d X %d.\n', n, n_rad, n_tan);                
                n=N+1;
            else
                n=n+1;
            end
end

%% comparison with Exact Solution

if RMS<RMS_acceptance
    for j=1:1:Jmax
        for i=1:1:Imax
            SUM=SUM+(PHI_Soln(i,j)-Exact_Soln(i,j))*(PHI_Soln(i,j)-Exact_Soln(i,j));
        end
    end

            RMS=sqrt(SUM/((Imax-2)*(Jmax-2)));
            fprintf('%d X %d ----->  RMS=%d  ----->  No. of iteration=%d\n', n_rad, n_tan, RMS, Iteration)
end


