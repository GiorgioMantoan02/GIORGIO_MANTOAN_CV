% Set dei parametri

clear all
clc
S_0=10;
K=9;
r=0.02;
T=2;
sigma=0.05;

%Discretizzazione
N=15;

%Calcolo con albero CRR
Call=european_option(S_0,sigma,r,K,T,N,1)
Put=european_option(S_0,sigma,r,K,T,N,-1)

%Calcolo analitico
[C,P]=blsprice(S_0,K,r,T,sigma)