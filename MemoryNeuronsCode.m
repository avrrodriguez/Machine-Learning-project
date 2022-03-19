

%%% Altered Spike-Time-Dependent Plasticity Following Impaired 
%%% Myelination of Thalamocortical Projections to Primary Motor Cortex

clear all

%Motor Cortex Excitatory Neurons 
NumM1E=800;
a(1:NumM1E,1)=0.02;
b(1:NumM1E,1)=0.2;
c(1:NumM1E,1)=-65+(15.*rand(NumM1E,1).^2);
d(1:NumM1E,1)=8-(6.*rand(NumM1E,1).^2);

%Motor Cortex Inhibitory Neurons
NumM1I=200;
All_M1nodes=NumM1E+ NumM1I;
a(NumM1E+1:All_M1nodes,1)=0.02+(0.08.*rand(NumM1I,1));
b(NumM1E+1:All_M1nodes,1)=0.2+(0.05.*rand(NumM1I,1));
c(NumM1E+1:All_M1nodes,1)=-65;
d(NumM1E+1:All_M1nodes,1)=2;

%Memory Cortex Excitatory Neurons
NumMemoryE=800;
a(All_M1nodes+1:All_M1nodes+NumMemoryE,1)=0.02;
b(All_M1nodes+1:All_M1nodes+NumMemoryE,1)=0.2;
c(All_M1nodes+1:All_M1nodes+NumMemoryE,1)=-65+(15.*rand(NumMemoryE,1).^2);
d(All_M1nodes+1:All_M1nodes+NumMemoryE,1)=8-(6.*rand(NumMemoryE,1).^2);

%Memory Cortex Inhibitory Neurons
NumMemoryI=200;
All_Memorynodes=NumMemoryE+ NumMemoryI;
a(All_M1nodes+NumMemoryE+1:All_M1nodes+All_Memorynodes,1)=0.02+(0.08.*rand(NumMemoryI,1));
b(All_M1nodes+NumMemoryE+1:All_M1nodes+All_Memorynodes,1)=0.2+(0.05.*rand(NumMemoryI,1));
c(All_M1nodes+NumMemoryE+1:All_M1nodes+All_Memorynodes,1)=-65;
d(All_M1nodes+NumMemoryE+1:All_M1nodes+All_Memorynodes,1)=2;

%Initialize weight matrix
All_nodes=All_M1nodes+All_Memorynodes;
SynWeights=rand(All_nodes);

%Column is the type
E_Scalar=0.5;%%
I_Scalar=-1;%%%

SynWeights(:,1:NumM1E)=SynWeights(:,1:NumM1E).*E_Scalar;%col 1-800 excitatory
SynWeights(:,All_M1nodes+1:All_M1nodes+NumMemoryE)=SynWeights(:,All_M1nodes+1:All_M1nodes+NumMemoryE).*E_Scalar;%col 1000-1800 excitatory
SynWeights(:,NumM1E+1:All_M1nodes)=SynWeights(:,NumM1E+1:All_M1nodes).*I_Scalar;%col 801-1000 inhibatory
SynWeights(:,All_M1nodes+NumMemoryE+1:All_nodes)=SynWeights(:,All_M1nodes+NumMemoryE+1:All_nodes).*I_Scalar;%col 1801-2000 inhibatory

%Deleting self-to-self weights
IDmat=eye(All_nodes);
IDRmat=(IDmat.*-1)+1;
SynWeights=IDRmat.*SynWeights;

%Deleting all cross network impacts
for row=1:All_nodes
    for col=1:All_nodes
        if (col>=All_M1nodes+1)&&(row<= All_M1nodes)
            SynWeights(row,col)=0;
        end
        if (col<=All_M1nodes && col>NumM1E)&&(row > All_M1nodes)
            SynWeights(row,col)=0;
        end
        if ((row >= All_M1nodes+NumMemoryE) && (row <=All_nodes)) && col <= All_M1nodes
            SynWeights(row,col) = 0;
        end
    end
end

%initialize VM and u as array
Vm=-65.*ones(All_nodes,1);
u=b.*Vm;

%Misc var.
Spikes=[];
dt=.5;
SimTime=1000;
PSP =zeros(All_nodes,1);
Vm_rec=[];
AP_Max=30;

%Strength of input
Ein=3.6;
Iin=1.8;


Spike_rec=zeros(All_nodes,SimTime);
for t=1:SimTime
    %rand PSP input to each node
    PSP=[Ein.*randn(NumM1E,1); Iin.*randn(NumM1I,1);2.9.*randn(NumMemoryE,1); 1.3.*randn(NumMemoryI,1)];
    
    SpikingNodes = find(Vm>=30);
    
    if(~isempty(SpikingNodes))
        Spike_rec(SpikingNodes,t)=1;
        Vm(SpikingNodes)=c(SpikingNodes);
        u(SpikingNodes)= u(SpikingNodes) + d(SpikingNodes);
        Vm_rec(SpikingNodes,t)=AP_Max;
        Spikes=[Spikes; t+(0.*SpikingNodes),SpikingNodes];
        PSP = PSP+ sum(SynWeights(:,SpikingNodes),2);
    end
    
    Vm = Vm+dt*((((0.04.*Vm.^2)+(5.*Vm)+140)-u)+PSP);
    Vm = Vm+dt*((((0.04.*Vm.^2)+(5.*Vm)+140)-u)+PSP);
    u = u+(a.*((b.*Vm)-u));
    Vm_rec(:,t)=Vm;
end

Simulations=6;

for SimLoop=1:Simulations
RandomMemoryENeuron = round((All_M1nodes+1) + ((All_M1nodes+NumMemoryE)-(All_M1nodes+1)).*rand(1));
PostN=Spike_rec(RandomMemoryENeuron,:);
PreN=[Spike_rec(All_M1nodes+1:RandomMemoryENeuron-1,:);Spike_rec(RandomMemoryENeuron+1:All_M1nodes+NumMemoryE,:)];
PreNSynWeights = rand(799,1) - 0.5;
[Pre_h,Pre_w]=size(PreN);

hLTtime=25;
PostSpikeTime=hLTtime+1;
tau_LTP=10;
tau_LTD=10.8;
Lrate=0.1;
t_stop=1000;

for PreNeuron=1:Pre_h
    SpikeCount=0;
    STDPw_rec=[];

        for t=hLTtime:t_stop-hLTtime
            if (PostN(1,t)==1)%PostN spike train
                SpikeCount= SpikeCount+1;
                PreSynWindow = PreN(PreNeuron, t-hLTtime:t+hLTtime); 
                PreSpikeTimes=find (PreSynWindow==1);
                tD=PostSpikeTime-PreSpikeTimes;
                if (isempty(tD))
                    STDPw=0;
                else
                    Mu_tD=mean(tD);

                    if (Mu_tD>=0)
                        STDPw=exp(-Mu_tD./tau_LTP);
                    else
                        STDPw=-exp(Mu_tD./tau_LTD).*0.5;
                    end
                end
                STDPw_rec(1,SpikeCount)=STDPw;
            end
        end
    DeltaW(SimLoop,PreNeuron)=mean(STDPw_rec);
    PreNSynWeights(PreNeuron,1)= PreNSynWeights(PreNeuron,1) + Lrate .* DeltaW(SimLoop,PreNeuron);

end
end

figure;
for SimLoop = 1:Simulations
subplot(2,3,SimLoop)
histogram(DeltaW(SimLoop,:),20);
xlabel('W');
ylabel('Frequency');
end

plot(Spikes(:,1),Spikes(:,2),'.');
