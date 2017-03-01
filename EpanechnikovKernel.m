function [ K ] = EpanechnikovKernel( dimension, x_metric )
%EpanechnikovKernel 根据输入的维数和x的度量值计算核函数值
K = ( (2.*pi).^(-dimension/2) ).*exp( -1/2.*(x_metric.^2) );

end

