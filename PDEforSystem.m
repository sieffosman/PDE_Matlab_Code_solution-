%%%%%%%To Solve our PDE
%u_t= (beta/4)*(s^2 -1)*u_ss-(1/eta)*(u_ssss)
close all 
clear all 

L=1;
N=50; %%%should be 51

%define eta and beta - add in constants 
beta=0.00617;
eta=1;

% Discretise over s
s= linspace(-L,L,N); 
ds=(s(2)-s(1)); %%% can think of as delta x 
% Discretise over t
dt=.000000;
t_max=1;
t=0:dt:t_max; 

%%%assign memory for displacement/pert
U =zeros(length(s),length(t));
size(U)

%%%Make derivatives look better:

U(:,1)=0.01*s.^2; %%small pertubation 


A=beta./4;
A=100;
B=1./eta;
B=1;
%loop over time 

for k = 1:length(t)-1 %%%integrates to k+1 - loop for time 
    %loop over space 
    for j = 3:length(s)-2 %%%%%%%%%%%% had to change to start @ 4 for BCs and end at length(s)-3
        
        %%%%To make Clearer let
        u_ss=U(j+1,k)-2*(U(j,k))+U(j-1,k);
        u_ss=u_ss./((ds).^2);
        
        u_ssss=U(j+2,k)-4*U(j+1,k)+6*U(j,k)-4*U(j-1,k)+U(j-2,k);
        u_ssss=u_ssss./((ds).^4);
        
        %%actual derivative expression 
        U(j,k+1)=U(j,k)+dt*(A.*((s(j).^2)-1).*(u_ss)-B.*(u_ssss));
        
     
        
        
        
       
        
       
    
    end
    
    
       %%Boundary conditions 
%%U_sss=U_ss=0 at s=+-1 for now use one from vid (Dr Carlos Mon...)
%U(space,time)=U(j,k)

%%Backwards BC for s=-1 (1st order from wiki formula)
%%U_ss(-1)=0 ....
U(2,k+1)=2*U(3,k+1)-U(4,k+1);
%%U_sss(-1)=0
U(1,k+1)=3*U(2,k+1)-3*U(3,k+1)+U(4,k+1);
%%Forwards BC for s=1 (1st order from wiki formula)
%%U_ss(1)=0
U(end-1,k+1)=2*U(end-2,k+1)-U(end-3,k+1);
%%U_sss(1)=0
U(end,k+1)=3*U(end-1,k+1)-3*U(end-2,k+1)+U(end-3,k+1);


              %%%end of BC stuff 
              
             
end

U;

%%%plot this 
hold on
axis equal 
for i=1:10000:length(U(1,:))
    plot(s,U(:,i)')
    pause 
   
end





%%This makes noise when code is done. 
load handel
sound(y,Fs)


%beep on;
%beep;

