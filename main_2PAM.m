%% noms des binomes

%%Initialisation:

clear ;
close all ;
clc

TEB = [];
teb=0;
%% Initialisation des variables
Fe=10000;
Te=1/10000;
Ds=1000;
Ts=1/1000;
Ns=5000;
Fse=Ts/Te;
N=512;
T_fin_axis=(50*Ts-Te)*Fe;
axe_teb_db = 0:0.5:10;
axe_teb = 10.^(axe_teb_db/10)
%% Emetteur

sb=randi([0,1],1,Ns);

ss=[];
for i=1:Ns
    if(sb(i)==0)
        ss=[ss -1];
    else
        ss=[ss 1];
    end
end

g=[1 1 1 1 1 1 1 1 1 1];

s=upsample(ss,10);
sl=conv(s,g);

for x=axe_teb_db;
    c=0;
    teb=0;
    
while (teb<100)
    %% Canal
    hl=[1];
    
    var = Fse/(2* 10^(x/10));
    nl=sqrt(var).*randn(1,length(sl)) ;
    yl=conv(sl,hl) +nl;
    
    %% R�cepteur
    ga=[1 1 1 1 1 1 1 1 1 1];
    
    rl=conv(ga,yl);
    
    rl_n=rl(Fse:Fse:Ns*10+1);
    
    
    S=pamdemod(rl_n,2);
    
    res=[];
    for i=1:Ns
        if(S(i)==0)
            res=[res 0];
        else
            res=[res 1];
        end
    end
    

    %% Calcul du TEB :
    nb=0;
    for i=1:length(sb)
        if (sb(i) ~= res(i))
            nb=nb+1;
        end
    end
    teb=teb+nb;
    c=c+1;
    
end
    teb = teb / (c*length(res));
    TEB = [TEB,teb];
    
end
   




%% figures


figure;
plot(sl);
xlabel("temps en sec")
ylabel("information en symbole")
xlim([0,T_fin_axis])
ylim([0,1])
grid on;

figure;
plot(rl);
xlabel("temps en sec")
ylabel("signal apres filtrage de reception")
xlim([0,T_fin_axis])
ylim([0,1])
grid on;

%% Calcul de la DSP th�orique :

axe_freq = 0:2500/N:2500;
DSP_sl = Ts * (sin(pi*axe_freq*Ts).^2) ./ (pi*axe_freq*Ts).^2 ;
DSP_ss=var/Ts*ones(1,2500);
figure;
semilogy(DSP_sl);
hold on;
semilogy(DSP_ss);
xlabel("frequence en Hz")
ylabel("DSP theoriques")
xlim([0,T_fin_axis])

%% Affichage de la DSP théorique :

figure;
pwelch(ss,N);
hold on;
pwelch(sl,N);

%% Affichage du diagramme de loeil

eyediagram(rl,30,3*Ts)


%% Evolution du TEB :


pb=0.5 * erfc(sqrt(axe_teb));

figure;
semilogy(10*log10(axe_teb),TEB);
hold on;
semilogy(10*log10(axe_teb),pb);
ylabel("Pb et TEB")
xlabel("Eb/No en db")
xlim([0 10]);





