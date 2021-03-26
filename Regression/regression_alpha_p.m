clc
clear all
close all

myDir        = './data';
myFiles      = dir(fullfile(myDir,'*.mat'));
alpha        = [];
date         = [];
dwelltime    = [];
dwell_time   = [];
length_dates = [];
number_dwell = [];
for k = 1:length(myFiles)
  baseFileName = myFiles(k).name;
  fullFileName = fullfile(myDir, baseFileName);
  if contains(baseFileName, 'alphas')
      load(strcat('./data/', baseFileName));
      alpha = cat(2, alpha, alphas);
      if all(alphas < 1)
         fprintf('Valid. All alpha smaller than 1. \n')
      else
         fprintf('Invalid. Exist alpha greater than 1. \n')
      end
  else 
      load(strcat('./data/', baseFileName));
      dates        = datenum(dates);
      length_dates = cat(1, length_dates, length(dates));
      dates_shift  = [dates(2:end); 0];
      dwell        = dates_shift - dates;
      dwell_time   = dwell(1:length(dwell)-1);
      number_dwell = cat(1, number_dwell, length(dwell_time));
      dwelltime    = cat(1, dwelltime, dwell_time);
  end 
end

%% find the index of the regions
test_data          = readtable('./dati-regioni/dpc-covid19-ita-regioni.csv');
region_index       = [];
region_index(:,1)  = find(contains(test_data.denominazione_regione,'Umbria'));
region_index(:,2)  = find(contains(test_data.denominazione_regione,'Marche'));
region_index(:,3)  = find(contains(test_data.denominazione_regione,'Lazio'));
region_index(:,4)  = find(contains(test_data.denominazione_regione,'Abruzzo'));
region_index(:,5)  = find(contains(test_data.denominazione_regione,'Molise'));
region_index(:,6)  = find(contains(test_data.denominazione_regione,'Campania'));
region_index(:,7)  = find(contains(test_data.denominazione_regione,'Puglia'));
region_index(:,8)  = find(contains(test_data.denominazione_regione,'Basilicata'));
region_index(:,9)  = find(contains(test_data.denominazione_regione,'Calabria'));
region_index(:,10) = find(contains(test_data.denominazione_regione,'Sicilia'));
region_index(:,11) = find(contains(test_data.denominazione_regione,'Piemonte'));
region_index(:,12) = find(contains(test_data.denominazione_regione,'Sardegna'));
region_index(:,13) = find(contains(test_data.denominazione_regione,"Valle d'Aosta"));
region_index(:,14) = find(contains(test_data.denominazione_regione,'Lombardia'));
region_index(:,15) = find(contains(test_data.denominazione_regione,'P.A. Bolzano'));
region_index(:,16) = find(contains(test_data.denominazione_regione,'Veneto'));
region_index(:,17) = find(contains(test_data.denominazione_regione,'Friuli Venezia Giulia'));
region_index(:,18) = find(contains(test_data.denominazione_regione,'Liguria'));
region_index(:,19) = find(contains(test_data.denominazione_regione,'Emilia-Romagna'));
region_index(:,20) = find(contains(test_data.denominazione_regione,'Toscana'));
region_index(:,21) = find(contains(test_data.denominazione_regione,'P.A. Trento'));

%% start index for each region
start_index(1)     = 83;  population(1)  = 880285;
start_index(2)     = 86;  population(2)  = 1518400;
start_index(3)     = 60;  population(3)  = 5865544;
start_index(4)     = 105; population(4)  = 1305770;
start_index(5)     = 76;  population(5)  = 302265;
start_index(6)     = 113; population(6)  = 5785861;
start_index(7)     = 87;  population(7)  = 4008296;
start_index(8)     = 77;  population(8)  = 556934;
start_index(9)     = 103; population(9)  = 1924701;
start_index(10)    = 77;  population(10) = 4968410;
start_index(11)    = 60;  population(11) = 4341375;
start_index(12)    = 102; population(12) = 1630474;
start_index(13)    = 59;  population(13) = 125501;
start_index(14)    = 59;  population(14) = 10103969;
start_index(15)    = 61;  population(15) = 1074819; %population of Trento[21] included
start_index(16)    = 68;  population(16) = 4907704;
start_index(17)    = 78;  population(17) = 1211357;
start_index(18)    = 87;  population(18) = 1543127;
start_index(19)    = 75;  population(19) = 4467118;
start_index(20)    = 64;  population(20) = 3722729;
start_index(21)    = 61;

for j = 1:20
    for i = 1:number_dwell(j)
        test_index(i,j) = region_index( start_index(j) ...
                        + sum( dwelltime( sum(number_dwell(1:j-1)) + 1 : sum(number_dwell(1:j-1)) + i)) ,j);
    end
end
test_index = cat(1, zeros(1,20) ,test_index);
for i = 1:20
   test_index(1,i) = region_index(start_index(i), i); 
end    
for i =  1:20
    people_tested(1:length_dates(i),i) = test_data.casi_testati(test_index(1:length_dates(i),i));
    people_tested_subtractor(1:length_dates(i),i) = test_data.casi_testati(test_index(1:length_dates(i),i)-21);
end
people_tested(1:length_dates(15),15) = people_tested(1:length_dates(15),15) + test_data.casi_testati(test_index(1:length_dates(15),15) + 1);%add the data from Bolzano and Trento
people_tested_subtractor(1:length_dates(15),15) = people_tested_subtractor(1:length_dates(15),15) + + test_data.casi_testati(test_index(1:length_dates(15),15) + 1);

for j = 1:20
    for i = 1:length_dates(j)-1
        people_tested_each_dwell(i,j) = people_tested(i+1,j) - people_tested_subtractor(i,j);
    end
end

for i = 1:20
    for j = 1:number_dwell(i)
        fraction_of_people_tested_daily( j+sum(number_dwell(1:i-1)) ) = people_tested_each_dwell(j,i) / dwelltime( sum(number_dwell(1:i-1))+j) / population(i);
    end
   
end 

alpha_matrix = zeros(50,20);
for i = 1:20
    for j = 1:number_dwell(i)
        alpha_matrix(j,i) = alpha( j+sum(number_dwell(1:i-1)) );
    end
end

%% regression in national level
ind = find(fraction_of_people_tested_daily<0);
fraction_of_people_tested_daily(ind) = [];
alpha(ind) = [];
regression_para = polyfit(fraction_of_people_tested_daily,alpha,1);
alpha_reg = polyval(regression_para,fraction_of_people_tested_daily);
figure
plot(fraction_of_people_tested_daily(:),alpha(:),'o')
hold on
plot(fraction_of_people_tested_daily(:),alpha_reg(:),'-')
xlabel('p')
ylabel('$\alpha$','interpreter','latex')
grid on

