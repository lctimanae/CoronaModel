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

classdef BergeronModel_ini
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
        function obj = BergeronModel_ini(values)
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

        function tao = taoLine(obj, vk_t_deltat , vm_t_deltat)
            v_t_deltat = (vk_t_deltat + vm_t_deltat) / 2;
            if v_t_deltat<obj.V_crit
                tao = sqrt(obj.L0*obj.C0)*obj.L;
            else
                C = obj.C0 + 2*obj.K_C*(1 - obj.V_crit/v_t_deltat);
                tao = sqrt(obj.L0*C)*obj.L;
            end
        end

        function Z = ZInitialLine(obj)
            R = obj.R0*obj.L;
            Z0 = sqrt(obj.L0/obj.C0);
            Z = Z0+R/2;
        end

        function [Ik_t_tao, Im_t_tao, Z]=historicalCurrents(obj , imk_t_tao , ikm_t_tao , vm_t_tao , vk_t_tao , vk_t_deltat , vm_t_deltat)
            v_t_deltat = (vk_t_deltat + vm_t_deltat) / 2;
            R = obj.R0*obj.L;

            if v_t_deltat<obj.V_crit
                Z0 = sqrt(obj.L0/obj.C0);
                Z = Z0+R/2;
                H=(Z0-R/2)/(Z0+R/2);
                Ik_t_tao = -(1/Z)*vm_t_tao - H*imk_t_tao;
                Im_t_tao = -(1/Z)*vk_t_tao - H*ikm_t_tao;
            else
                C = obj.C0 + 2*obj.K_C*(1 - obj.V_crit/v_t_deltat);
                Z_0k = sqrt(obj.L0/C);
                G = obj.K_G * (1 - obj.V_crit/v_t_deltat)^2;
                Rs = 1/(G*obj.L);
                %Source: T. M. Pereira and M. C. Tavares, "Development of a Voltage-Dependent Line Model to Represent the Corona Effect in Electromagnetic Transient Program," in IEEE Transactions on Power Delivery, vol. 36, no. 2, pp. 731-739, April 2021
                Temp = 2*Rs*R + 4*Rs*Z_0k + R*Z_0k;
                Z = Temp/(4*Rs + 2*Z_0k);
                % Corrected equations, proposed by (Pereira and Tavares, 2021)
                Ik_t_tao =  ((2*Z_0k-4*Rs)/Temp)*vm_t_tao +((4*Rs*R)/Temp - 1)*imk_t_tao;
                Im_t_tao =  ((2*Z_0k-4*Rs)/Temp)*vk_t_tao +((4*Rs*R)/Temp - 1)*ikm_t_tao;
                % Equivalent, non-simplified formulas
                %Ik_t_tao = -(4*Rs*(vm_t_tao - (R*imk_t_tao)/2 + Z_0k*(imk_t_tao - (vm_t_tao - (R*imk_t_tao)/2)/(2*Rs))))/(2*Rs*R + 4*Rs*Z_0k + R*Z_0k);
                %Im_t_tao = -(4*Rs*(vk_t_tao - (R*ikm_t_tao)/2 + Z_0k*(ikm_t_tao - (vk_t_tao - (R*ikm_t_tao)/2)/(2*Rs))))/(2*Rs*R + 4*Rs*Z_0k + R*Z_0k);
            end
        end

    end
end