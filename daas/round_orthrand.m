function [x]=round_orthrand(tt,varargin)
  %Approximate TT-tensor with another one with specified accuracy
  %   [TT]=ROUND(TT,EPS) Approximate TT-tensor with relative accuracy EPS
  %
  %   [TT]=ROUND(TT,EPS,l) Approximate TT-tensor with relative accuracy 
  %   EPS and maximal rank l. l can be array of ranks or a number
  %
  %
  % TT-Toolbox 2.2, 2009-2012
  %
  %This is TT Toolbox, written by Ivan Oseledets et al.
  %Institute of Numerical Mathematics, Moscow, Russia
  %webpage: http://spring.inm.ras.ru/osel
  %
  %For all questions, bugs and suggestions please mail
  %ivan.oseledets@gmail.com
  %---------------------------
  
  %
  if nargin == 2
    l = varargin{1};
  elseif nargin == 1
    l = zeros(d + 1, 1);
    l(1) = 1; l(end) = 1;
    for i = 1 : d - 1
      L = prod(n(1 : i));
      R = prod(n(i + 1 : end));
      l(i + 1) = min(L, R);
    end
  else   
    error("Invalid number of arguments.");
  end
  d=tt.d;
  n=tt.n;
  r=tt.r;
  if (isscalar(l) && length(n) > 1)
    l = [1, ones(1, length(n)-1)*l, 1]';
  end
  if length(l) ~= 1 && (l(1) ~= 1 || l(end) ~= 1)
    error('The first and the last TT-ranks must be equal to 1');  
  else
    for i = 1 : length(n) - 1
      L = prod(n(1 : i));
      R = prod(n(i + 1 : end));
      l(i + 1) = min([L, R, l(i + 1)]);
    end    
  end
  pos=tt.ps;
  cr=tt.core;
  pos1=1;
  % record the norm of each R of QR decomposition
  nrm=zeros(d,1);
  core0=cr(1:r(1)*n(1)*r(2));
  %Orthogonalization from left-to-tight
  for i=1:d-1
    % V reshape left core
    core0=reshape(core0,[r(i)*n(i),r(i+1)]);
    % Organize the left core
     [core0,ru]=qr(core0,0); 
     nrm(i+1)=norm(ru,'fro');
     if (nrm(i+1)~=0)
      ru=ru./nrm(i+1);
     end;
    % get the right core and left-multiply it with ru
    core1=cr(pos(i+1):pos(i+2)-1);
    core1=reshape(core1,[r(i+1),n(i+1)*r(i+2)]);
    core1=ru*core1;
    r(i+1)=size(core0,2);
    cr(pos1:pos1-1+r(i)*n(i)*r(i+1))=core0(:);
    cr(pos1+r(i)*n(i)*r(i+1):pos1+r(i)*n(i)*r(i+1)+r(i+1)*n(i+1)*r(i+2)-1)=core1(:);
    core0=core1;
    pos1=pos1+r(i)*n(i)*r(i+1);
  end
  pos1=pos1+r(d)*n(d)*r(d+1)-1;
  cr=cr(1:pos1); %Truncate storage if required
  pos=cumsum([1;n.*r(1:d).*r(2:d+1)]); 
  core0=cr(pos1-r(d)*n(d)*r(d+1)+1:pos1);
  x = tt_tensor;
  ps_x = cumsum([1; n(1 : end) .* l(1 : end - 1) .* l(2 : end)]);
  cr_x = zeros(ps_x(end) - 1, 1);
  pos_x = ps_x(end);
  % r_x = ones(d + 1, 1);
  for i=d:-1:2
    core1=cr(pos(i-1):pos(i)-1); 
    core0=reshape(core0,r(i),n(i) * l(i + 1));
    core1=reshape(core1,r(i-1)*n(i-1),r(i));

    R = randn(l(i), r(i));
    Y = R * core0;
    [Q, ~] = qr(Y', 0);
    core1 = core1 * core0 * Q;
    core0 = Q';
    l(i) = size(core0, 1);
    % try
    %   [u,s,v]=svd(core0, 'econ');
    % catch err
    %   if (strcmp(err.identifier,'MATLAB:svd:svdNoConvergence'))
    %     % When svd doesn't converge, svds usually works fine (but much slower).
    %     warning('SVD did not converge, switched to svds.');
    %     [u,s,v]=svds(core0, min(size(core0)));
    %   else
    %     rethrow(err);
    %   end
    % end

    % s=diag(s);
    % r1=my_chop2(s,norm(s)*ep);
    % r1=min(r1,l(i));

    % u=u(:,1:r1);s=s(1:r1); v=v(:,1:r1);
    % u=u*diag(s);
    % r(i)=r1;
    % core1=core1*u;
    % core0=v';
    cr_x(pos_x - l(i) * n(i) * l(i + 1) : pos_x - 1)=core0(:);
    pos_x = pos_x - l(i) * n(i) * l(i + 1);
    % r_x(i) = size(core0, 1);
    
    % cr(pos1-r(i)*n(i)*r(i+1)-r(i-1)*n(i-1)*r(i)+1:pos1-r(i)*n(i)*r(i+1))=core1(:);
    % pos1=pos1-r(i)*n(i)*r(i+1);
    core0=core1;
  end
  cr_x(pos_x-l(1)*n(1)*l(2):pos_x - 1)=core0(:);
  pos_x = pos_x - l(1)*n(1)*l(2);
  cr_x = cr_x(pos_x:end); %Truncate unwanted elements;
  x.r=l;
  % x.core = cr_x;
  x.n = n;
  x.d = d;
  x.ps = cumsum([1; n .* l(1 : end - 1) .* l(2 : end)]);
  pp=cr_x(1:x.ps(2)-1);
  nrm(1) = norm(pp,'fro');
  if (nrm(1) ~= 0)
      pp = pp./nrm(1);
  end;
  cr_x(1:l(1)*n(1)*l(2)) = pp;
   %Now a simple trick: balance the product of numbers;
   %All cores are orthogonal except the first one. Thus, we know the norm
   nrm0=sum(log(abs(nrm))); 
   nrm0=nrm0/d; nrm0=exp(nrm0);
   ps=x.ps;
   for i=1:d
      cr_x(ps(i):ps(i+1)-1)=nrm0*cr_x(ps(i):ps(i+1)-1);
   end
   x.core=cr_x;
   % tt.over=0;
  return
  end
  