function [ output_args ] = graficoIpIfDezDep( temp , shot)

figure('Name',strcat(['Sparo ',num2str(shot)]))
subplot(3,1,1)
plot(temp.mytime,temp.IPLmis, 'LineWidth',2); grid on;
grid on; legend('Ip');
xlabel('tempo[s]','interpreter','latex','fontsize',12);
ylabel('[$A$]','interpreter','latex','fontsize',12);
title('Plasma current.');
subplot(3,1,2)
plot(temp.mytime,temp.Ifm, 'LineWidth',2); grid on;
grid on; legend('If');
xlabel('tempo[s]','interpreter','latex','fontsize',12);
ylabel('[$A$]','interpreter','latex','fontsize',12);
title('Plasma current.');
subplot(3,1,3)
plot(temp.mytime,temp.dIfm, 'LineWidth',2); grid on;
grid on; legend('DIf');
xlabel('tempo[s]','interpreter','latex','fontsize',12);
ylabel('[$\frac{A}{s}$]','interpreter','latex','fontsize',12);
title('dIfm');

figure('Name',strcat(['DezDepRs2Rs1 Sparo ',num2str(shot)]))
subplot(2,1,1)
plot(temp.mytime,temp.dep,temp.mytime,temp.dez, 'LineWidth',2); grid on;
grid on; legend('dep','dez');
xlabel('tempo[s]','interpreter','latex','fontsize',12);
ylabel('Web','interpreter','latex','fontsize',12);
title('Flusso magnetico');
subplot(2,1,2)
plot(temp.mytime,temp.rs2,...
     temp.mytime,temp.myRS2prep,...
     temp.mytime,temp.rs1,...
     temp.mytime,temp.myRS1prep, 'LineWidth',2); grid on;
grid on; legend('rs2','rs2 prep','rs1','rs1 prep');
xlabel('tempo[s]','interpreter','latex','fontsize',12);
ylabel('Web','interpreter','latex','fontsize',12);
title('rs2/rs1');



end

