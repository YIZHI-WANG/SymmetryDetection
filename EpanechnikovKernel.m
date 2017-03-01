function [ K ] = EpanechnikovKernel( dimension, x_metric )
%EpanechnikovKernel ���������ά����x�Ķ���ֵ����˺���ֵ
K = ( (2.*pi).^(-dimension/2) ).*exp( -1/2.*(x_metric.^2) );

end

