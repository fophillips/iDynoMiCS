package simulator.reaction.kinetic;

import org.jdom.Element;
import utils.LogFile;
import utils.XMLParser;

/**
 * Created with IntelliJ IDEA.
 * Date: 01/07/2014
 * Time: 11:42
 *
 * @author Fred Phillips
 */
public class AcidKinetic extends IsKineticFactor
{
	/**
	 * Temperature
	 */
	private Double _T;
	/**
	* Rate constant for proton-promoted dissolution
	*/
	private Double _kH;
	/**
	 * Reaction order of proton-promoted dissolution
	 */
	private Double _nH;
	/**
	* Rate constant for hydroxide-promoted dissolution
	*/
	private Double _kOH;
	/**
	 * Reaction order of proton-promoted dissolution
	 */
	private Double _nOH;
	/**
	 * Activation energy
	 */
	private Double _Ea;
	/**
	 * Molar mass
	 */
	private Double _molMass;
	/**
	 * Gas constant
	 */
	private static Double _R = 8.314; // J K−1 mol−1
	/**
	 * The standard amount concentration for calculating the activity
	 */
	private static Double _standardConc = 1.0; // mol L-1

	@Override
	public void init(Element defMarkUp)
	{
		_T = (new XMLParser(defMarkUp)).getParamDbl("T");
		_kH = (new XMLParser(defMarkUp)).getParamDbl("kH");
		_nH = (new XMLParser(defMarkUp)).getParamDbl("nH");
		_kOH = (new XMLParser(defMarkUp)).getParamDbl("kOH");
		_nOH = (new XMLParser(defMarkUp)).getParamDbl("nOH");
		_Ea = (new XMLParser(defMarkUp)).getParamDbl("Ea");
		_molMass = (new XMLParser(defMarkUp)).getParamDbl("molMass");
		nParam = 7;
	}

	@Override
	public void initFromAgent(Element defMarkUp, double[] kineticParam, int paramIndex)
	{
		kineticParam[paramIndex] = (new XMLParser(defMarkUp)).getParamDbl("T");
		kineticParam[paramIndex+1] = (new XMLParser(defMarkUp)).getParamDbl("kH");
		kineticParam[paramIndex+2] = (new XMLParser(defMarkUp)).getParamDbl("nH");
		kineticParam[paramIndex+3] = (new XMLParser(defMarkUp)).getParamDbl("kOH");
		kineticParam[paramIndex+4] = (new XMLParser(defMarkUp)).getParamDbl("nOH");
		kineticParam[paramIndex+5] = (new XMLParser(defMarkUp)).getParamDbl("Ea");
		kineticParam[paramIndex+6] = (new XMLParser(defMarkUp)).getParamDbl("molMass");
		nParam = 7;
	}

	@Override
	public double kineticValue(double solute)
	{
		Double r = rate(solute, _kH, _nH, _kOH, _nOH, _Ea, _molMass, _T);
		return r;
	}

	@Override
	public double kineticDiff(double solute)
	{
		Double dR = diffRate(solute, _kH, _nH, _kOH, _nOH, _Ea, _molMass, _T);
		return dR;
	}

	@Override
	public double kineticValue(double solute, double[] paramTable, int index)
	{
		Double T = paramTable[index];
		Double kH = paramTable[index+1]; 
		Double nH = paramTable[index+2];
		Double kOH = paramTable[index+3];
		Double nOH = paramTable[index+4];
		Double Ea = paramTable[index+5];
		Double molMass = paramTable[index+6];

		Double r = rate(solute, kH, nH, kOH, nOH, Ea, molMass, T);

		return r;
	}

	@Override
	public double kineticDiff(double solute, double[] paramTable, int index)
	{
		Double T = paramTable[index];
		Double kH = paramTable[index+1]; 
		Double nH = paramTable[index+2];
		Double kOH = paramTable[index+3];
		Double nOH = paramTable[index+4];
		Double Ea = paramTable[index+5];
		Double molMass = paramTable[index+6];

		Double dR = diffRate(solute, kH, nH, kOH, nOH, Ea, molMass, T);

		return dR;
	}

	private Double rate(Double solute, Double kH, Double nH, Double kOH, Double nOH, Double Ea, Double molMass, Double T)
	{
		Double hPlus = hPlusConc(solute, molMass);
		Double rate;
		if(solute > 0)
		{
			rate = (kH * Math.pow(hPlus, nH) + kOH * Math.pow(hPlus, -nOH) * Math.pow(10, -14 * nOH))
					* Math.exp(-Ea/(_R * T));
 		}
		else
		{
			rate = 0.0;
		}
		return rate;
	}

	private Double diffRate(Double solute, Double kH, Double nH, Double kOH, Double nOH, Double Ea, Double molMass, Double T)
	{
		Double hPlus = hPlusConc(solute, molMass);
		Double dR;
		if(solute > 0)
		{
			dR = (kH * nH * Math.pow(hPlus, nH-1)
					- kOH * nOH * Math.pow(hPlus, -nOH - 1) * Math.pow(10,-14 * nOH))
					* Math.exp(-Ea/(_R * T));
		}
		else
		{
			dR = 0.0;
		}
		return dR;
	}

	private Double hPlusConc(Double solute, Double molMass)
	{
		return solute / (molMass * _standardConc);
	}
}

