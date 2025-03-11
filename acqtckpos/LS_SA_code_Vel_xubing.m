function [VR, dtRV, pr_rate] = LS_SA_code_Vel_xubing(XR_approx, XS, VS, dop_R, lambda, dtSV)

% SYNTAX:
%   [XR, dtR, cov_XR, var_dtR, PDOP, HDOP, VDOP, cond_num, bad_obs, bad_epoch, sigma02_hat, residuals_obs, is_bias] = LS_SA_code(XR_approx, XS, pr_R, snr_R, elR, distR_approx, dtS, err_tropo_RS, err_iono_RS, sys, SPP_threshold);
%
% INPUT:
%   XR_approx     = receiver approximate position (X,Y,Z)
%   XS            = satellite position (X,Y,Z)
%   pr_R          = code observations
%   snr_R         = signal-to-noise ratio
%   elR           = satellite elevation (vector)
%   distR_approx  = approximate receiver-satellite distance (vector)
%   dtS           = satellite clock error (vector)
%   err_tropo_RS  = tropospheric error
%   err_iono_RS   = ionospheric error
%   sys           = array with different values for different systems
%   SPP_threshold = maximum RMS of code single point positioning to accept current epoch
%
% OUTPUT:
%   XR   = estimated position (X,Y,Z)
%   dtR  = receiver clock error (scalar)
%   cov_XR  = estimated position error covariance matrix
%   var_dtR = estimated clock error variance
%   PDOP = position dilution of precision
%   HDOP = horizontal dilution of precision
%   VDOP = vertical dilution of precision
%   cond_num = condition number on the eigenvalues of the N matrix
%   bad_obs = vector with ids of observations found as outlier
%   bad_epoch = 0 if epoch is ok, -1 if there is no redoundancy, +1 if a posteriori sigma is greater than SPP_threshold
%   sigma02_hat = [a posteriori sigma (SPP sigma), v_hat'*(invQ)*v_hat), n-m] 
%   residuals_obs = vector with residuals of all input observation, computed from the final estimates
%   is_bias = inter-systems bias (vector with all possibile systems)
% DESCRIPTION:
%   Absolute positioning by means of least squares adjustment on code
%   observations. Epoch-by-epoch solution.
%----------------------------------------------------------------------------------------------
%                           goGPS v0.4.3
%
% Copyright (C) 2009-2014 Mirko Reguzzoni, Eugenio Realini
%
% Portions of code contributed by Stefano Caldera
%----------------------------------------------------------------------------------------------
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.
%----------------------------------------------------------------------------------------------

v_light = 299792458;

%number of observations
n = length(dop_R);

%number of unknown parameters
m = 4;

%approximate receiver-satellite distance
XR_mat = XR_approx(:,ones(n,1))'; 
distR_approx = sqrt(sum((XS-XR_mat).^2 ,2)); % a column vector

%design matrix
A = [(XR_approx(1) - XS(:,1)) ./ distR_approx, ... %column for X coordinate
     (XR_approx(2) - XS(:,2)) ./ distR_approx, ... %column for Y coordinate
     (XR_approx(3) - XS(:,3)) ./ distR_approx, ... %column for Z coordinate
      ones(n,1)];        %column for receiver clock delay (multiplied by c)
  
%known term vector  
b = sum(A(:,1:3).*VS,2) - dtSV' ; % m/s

%observation vector
y0 = -dop_R .* lambda; % m/s

%observation covariance matrix
% Q = cofactor_matrix_SA(elR, snr_R);
Q = eye(n);
invQ=diag((diag(Q).^-1));

%normal matrix
N = (A'*(invQ)*A);

%least squares solution
x = (N^-1)*A'*(invQ)*(y0-b);

VR  = x(1:3);
dtRV = x(4) ; % in unit of m/s

pr_rate = y0-b - dtRV ;


