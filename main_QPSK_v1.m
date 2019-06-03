%%Initialisation:

clear ;
close all ;
clc

TEB = [];

%% Initialisation des variables
Fe=10000;
Te=1/10000;
Ds=1000;
Ts=1/1000;
M=4;
fo=2500;
Ns=5000;
Fse=Ts/Te;
N=512;
axe_fourier=-Fe/2:Fe/2-Fe/N;
T_fin_axis=(50*Ts-Te)*Fe;
axe_teb_db = 0:0.5:10;
axe_teb = 10.^(axe_teb_db/10);

%% filtre racine cos sur-elevé

y=rcosdesign(0.5,8,Fse,'sqrt');

%% chaine de communication
%% Emetteur
for x=axe_teb_db
    var = Fse/(2*10^(x/10));
    
 %while (teb<100)

sb=randi([0,1],Ns/2,2);
sb_b=bi2de(sb);
ss=pskmod(sb_b,M,pi/M);

s=upsample(ss,10);
sl1=conv(s,y);    
axe_temp=(1:1:length(sl1))*Te;
sli=real(sl1)'.*cos(2*pi*fo*axe_temp);
slq=imag(sl1)'.*sin(2*pi*fo*axe_temp);
sl=sli+slq;

    
    
%% Canal bruité et echantillonage
axe_temp_recep=(1:1:length(sl))*Te;
nl=sqrt(var).*randn(1,length(sl)) ;
yl2 = zeros(1,length(sl));
for i=1:length(sl)
    yl2(i)= sl(i) + nl(i);
end
si=real(yl2).*cos(2*pi*fo*axe_temp_recep);
sq=-real(yl2).*sin(2*pi*fo*axe_temp_recep)*1i;
w=si+sq;

rl=conv(w,y);


rl_n=rl(length(y):Fse:end);
%% demodulation
S=pskdemod(rl_n,M,pi/M);

S_b=de2bi(S);

for l=1:length(S_b)
    if S_b(l,2)==S_b(l,1)
        if S_b(l,2)==0
            S_b(l,2)=1;
            S_b(l,1)=1;
        else
            S_b(l,2)=0;
            S_b(l,1)=0;
        end  
    else
    o=S_b(l,2);
    S_b(l,2)= S_b(l,1);
    S_b(l,1)=o; 
    end
end
        
%% Calcul du TEB :
nb=0;
for i=1:length(sb)
    if ((sb(i,1) ~= S_b(i,1)) || (sb(i,2) ~= S_b(i,2)))
        nb=nb+1;
        
    end
    
end

%   c=c+1
%  end
teb = nb / (length(sb));
TEB = [TEB,teb] ;   

end

pb=2 * qfunc(sqrt(2*log(M)*sin(pi/M)^2*axe_teb));

%% Figures

figure
plot(y);
xlabel("temps en sec")
ylabel("Reponse imp d'un filtre en racine de cosinus sur-élevé")

fvtool(y);

eyediagram(conv(s,y),3*Fse,3*Ts)
eyediagram(rl,3*Fse,3*Ts)
figure;
pwelch(sl,N);
% theor_sl=4*0.5*((cos((1.5*pi*axe_teb_db)./Ts)+sin((0.5*pi*axe_teb_db)./Ts)/(4*0.5*axe_teb_db./Ts))/(pi*sqrt(Ts*(1-(2*axe_teb_db./Ts).^2))));
% fft_sl=fft(theor_sl);
figure;
plot(real(rl));
xlim([0;50*Ts*Fe])
xlabel("temps en sec")
ylabel("partie réelle du signal modulé apres le filtre de reception")

figure;
pwelch(si,N);

figure;
semilogy(10*log10(axe_teb),TEB);
hold on;
semilogy(10*log10(axe_teb),pb);
xlim([0 10]);
ylabel("Pb et TEB")
xlabel("Eb/No en db")


