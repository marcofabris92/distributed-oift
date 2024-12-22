function [sg, dsig, ddsig] = sig(s,d,k)
%function [sg, dsig, ddsig] = sig(s,d,k)
%
%     s can be a vector of nonnegative values
%       that should be appropriately restricted depending
%       on which function is desired and what derivs are needed
%
%     d is a positive scalar
%
%     k = 0 or 1
%
%   break the function
%
%     sig(s)  =  ( sqrt( | s - d^2 | + d^2 ) - d )^2
%
%   into two functions (k=0,1)
%
%     sig0(s) = ( sqrt( 2*d^2 - s ) - d )^2     s < d^2
%
%     sig1(s) = ( sqrt( s )         - d )^2     s > d^2
%
%   which we hope can be glued together to give a C^2 function on s > 0 .
%
%   note that both of these functions are well defined for  0 <= s <= 2*d^2
%     and smooth (even analytic) on  0 < s < 2*d^2
%
%   thus we can compare the two functions (+ derivs) on that domain
%
%
% JH jul18 boulder


scale = 1/(d*(sqrt(2)-1))^2;
% scale

if 0 == k
  sqrt2d2ms = sqrt( 2*d^2 - s );

  % sg = ( sqrt( 2*d^2 - s ) - d ).^2;  % sig0
  sg   = (    sqrt2d2ms      - d ).^2;
  sg = scale*sg;

  if nargout > 1
    % dsig = d./sqrt(2*d^2 - s) - 1;    % dsig0
    dsig   = d./sqrt2d2ms       - 1;
    dsig = scale*dsig;

    if nargout > 2
      % ddsig = d./(2*(2*d^2 - s).^(3/2));    % ddsig0
      ddsig   = (d/2)./sqrt2d2ms.^3;
      ddsig = scale*ddsig;
    end

  end
elseif 1 == k
  r = sqrt(s);

  % sig = ( sqrt( s ) - d ).^2;    % sig1
  sg    = (    r      - d ).^2;
  sg = scale*sg;

  if nargout > 1
    % dsig = 1 - d./sqrt(s);    % dsig1   
    dsig   = 1 - d./r;
    dsig = scale*dsig;

    if nargout > 2
      % ddsig = d./(2*s.^(3/2));    % ddsig1
      ddsig   = (d/2)./r.^3;
      ddsig = scale*ddsig;
    end

  end
else
  error('k = %d is not 0 or 1',k);
end
