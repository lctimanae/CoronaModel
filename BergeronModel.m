%{ 
Author: Luis C. Timaná E.

R0: series resistance (Ohm/m)
L0: series inductance (H/m)
C0: shunt capacitance (F/m)
L: line length (m)
tao: travel time (s)
Z0: Characteristic impedance (Ohms)
Z: Norton interface impedance (Ohms)
R: total line resistance (Ohms)
Ik_t_tao, Im_t_tao: historical currents (A)
V_crit: Corona inception voltage (V) (value to be optimized)
K_C: Constant associated with shunt capacitance (value to be optimized)
K_G: Constant associated with shunt resistance (value to be optimized)

Reference:
- https://www.pscad.com/webhelp-v502-ol/EMTDC/Transmission_Lines/The_Bergeron_Model.htm
%}

classdef BergeronModel
    properties
        R0 {mustBeNumeric}
        L0 {mustBeNumeric}
        C0 {mustBeNumeric}
        L {mustBeNumeric}
        V_crit {mustBeNumeric}
        K_C {mustBeNumeric}
        K_G {mustBeNumeric}
    end

    methods
        % Constructor
        function obj = BergeronModel(values)
            if nargin == 1
                obj.R0 = values(1);
                obj.L0 = values(2);
                obj.C0 = values(3);
                obj.L = values(4);
                obj.V_crit = values(5);
                obj.K_C = values(6);
                obj.K_G = values(7);
            end
        end

        function C = CalculateCk(obj, vk)
            if vk<obj.V_crit
                C = obj.C0;
            else
                C = obj.C0 + 2*obj.K_C*(1 - obj.V_crit/vk);
            end
        end

        function tao = taoLine(obj, vk_t_deltat)
            C = CalculateCk(obj, vk_t_deltat);
            tao = sqrt(obj.L0*C)*obj.L;
        end

        function tao = taoLineSection(obj, vk_t_deltat , vkp1_t_deltat)
            % Method 1
            % C_k = CalculateCk(obj, vk_t_deltat);
            % C_kp1 = CalculateCk(obj, vkp1_t_deltat);
            % C = (C_k + C_kp1)/2;
            % tao = sqrt(obj.L0*C)*obj.L;

            % Method 2
            % tao_k = taoLine(obj, vk_t_deltat);
            % tao_kp1 = taoLine(obj, vkp1_t_deltat);
            % tao = (tao_k+tao_kp1)/2;

            % Method 3
            v_t_deltat = (vk_t_deltat + vkp1_t_deltat)/2;
            tao = taoLine(obj, v_t_deltat);
        end

        function Z = ZInitialLine(obj)
            R = obj.R0*obj.L;
            Z0 = sqrt(obj.L0/obj.C0);
            Z = Z0+R/2;
        end

        function Z = CalculateZLine(obj, vk_t_deltat)
            if vk_t_deltat<obj.V_crit 
                Z = ZInitialLine(obj);
            else
                R = obj.R0*obj.L;
                Z_0k = Calculate_Z(obj, vk_t_deltat);
                Rsk = Calculate_Rs(obj, vk_t_deltat);
                Z = (4*Rsk*Z_0k + 2*Rsk*R + R*Z_0k) / (4*Rsk + 2*Z_0k);
            end
        end

        function Z_0k=Calculate_Z(obj, v_k)
            if v_k<obj.V_crit
                Z_0k = sqrt(obj.L0/obj.C0);
            else
                C_k = CalculateCk(obj, v_k);
                Z_0k = sqrt(obj.L0/C_k);
            end
        end

        function Rs_k=Calculate_Rs(obj, v_k) % This result is valid only if v_k>=V_crit
            Gs_k = obj.K_G * (1 - obj.V_crit/v_k)^2;
            Rs_k = 1/(Gs_k*obj.L);
        end

        function [Ik, Im]=historicalCurrents(obj , ikm_t_taom , imk_t_taok , vk_t_taom , vm_t_taok  , vk_t_deltat , vm_t_deltat , vk_t_taom_deltat , vm_t_taok_deltat)
            R = obj.R0*obj.L;
            Z0 = sqrt(obj.L0/obj.C0);
            Z_0m_t_taok = Calculate_Z(obj, vm_t_taok_deltat);
            Z_0k_t_taom = Calculate_Z(obj, vk_t_taom_deltat);
            Rsm_t_taok = Calculate_Rs(obj, vm_t_taok_deltat);% This result is valid only if vm_t_taok_deltat>=obj.V_crit
            Rsk_t_taom = Calculate_Rs(obj, vk_t_taom_deltat);% This result is valid only if vk_t_taom_deltat>=obj.V_crit
            Rsk_t = Calculate_Rs(obj, vk_t_deltat);% This result is valid only if vk_t_deltat>=obj.V_crit
            Rsm_t = Calculate_Rs(obj, vm_t_deltat);% This result is valid only if vm_t_deltat>=obj.V_crit
            Z_0k = Calculate_Z(obj, vk_t_deltat);
            Z_0m = Calculate_Z(obj, vm_t_deltat);

            % Calculation of historical currents Ik
            if ( vk_t_deltat<obj.V_crit && vm_t_taok_deltat<obj.V_crit )
                Ik = vm_t_taok * ( -2 / (2*Z0+R)) + imk_t_taok * ( (R-2*Z0) / (R+2*Z0));
            elseif ( vk_t_deltat<obj.V_crit && vm_t_taok_deltat>=obj.V_crit )
                Ik = ( 2 / (Rsm_t_taok*(4*Z0+2*R)) ) * ( vm_t_taok*(-2*Rsm_t_taok + Z_0m_t_taok)  +  imk_t_taok*(R*Rsm_t_taok-2*Z_0m_t_taok*Rsm_t_taok-(R/2)*Z_0m_t_taok) );
            elseif ( vk_t_deltat>=obj.V_crit && vm_t_taok_deltat<obj.V_crit )
                Ik = ( (2*Rsk_t) / (4*Rsk_t*Z_0k+R*Z_0k+2*R*Rsk_t) ) * ( -2*vm_t_taok  +  imk_t_taok*(R-2*Z0)  );
            else
                temp = (2*Rsk_t) / ( Rsm_t_taok *(4*Rsk_t*Z_0k + R*Z_0k + 2*R*Rsk_t));
                Ik = temp * ( vm_t_taok*( -2*Rsm_t_taok + Z_0m_t_taok ) +  imk_t_taok * (R*Rsm_t_taok - 2*Z_0m_t_taok*Rsm_t_taok - (R/2)*Z_0m_t_taok));
            end

            % Calculation of historical currents Im
            if ( vm_t_deltat<obj.V_crit && vk_t_taom_deltat<obj.V_crit )
                Im = vk_t_taom * ( -2 / (2*Z0+R)) + ikm_t_taom * ( (R-2*Z0) / (R+2*Z0));
            elseif ( vm_t_deltat<obj.V_crit && vk_t_taom_deltat>=obj.V_crit )
                Im = ( 2 / (Rsk_t_taom*(4*Z0+2*R)) ) * ( vk_t_taom*(-2*Rsk_t_taom + Z_0k_t_taom)  +  ikm_t_taom*(R*Rsk_t_taom-2*Z_0k_t_taom*Rsk_t_taom-(R/2)*Z_0k_t_taom) );
            elseif ( vm_t_deltat>=obj.V_crit && vk_t_taom_deltat<obj.V_crit )
                Im = ( (2*Rsm_t) / (4*Rsm_t*Z_0m+R*Z_0m+2*R*Rsm_t) ) * ( -2*vk_t_taom  +  ikm_t_taom*(R-2*Z0)  );
            else
                temp = (2*Rsm_t) / ( Rsk_t_taom *(4*Rsm_t*Z_0m + R*Z_0m + 2*R*Rsm_t));
                Im = temp * ( vk_t_taom*( -2*Rsk_t_taom + Z_0k_t_taom ) +  ikm_t_taom * (R*Rsk_t_taom - 2*Z_0k_t_taom*Rsk_t_taom - (R/2)*Z_0k_t_taom));
            end

        end

    end
end