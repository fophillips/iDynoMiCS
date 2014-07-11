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
	 * Activation energy
	 */
	private Double _Ea;
	/**
	 * Acid molar mass
	 */
	private Double _acidMolarMass;
	/**
	 * Iron molar mass
	 */
	private Double _ironMolarMass;
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
		_T  = (new XMLParser(defMarkUp)).getParamDbl("T");
		_kH = (new XMLParser(defMarkUp)).getParamDbl("kH");
		_nH = (new XMLParser(defMarkUp)).getParamDbl("nH");
		_Ea = (new XMLParser(defMarkUp)).getParamDbl("Ea");
		_acidMolarMass = (new XMLParser(defMarkUp)).getParamDbl("acidMolarMass");
		_ironMolarMass = (new XMLParser(defMarkUp)).getParamDbl("ironMolarMass");
		nParam = 6;
	}

	@Override
	public void initFromAgent(Element defMarkUp, double[] kineticParam, int paramIndex)
	{
		kineticParam[paramIndex]   = (new XMLParser(defMarkUp)).getParamDbl("T");
		kineticParam[paramIndex+1] = (new XMLParser(defMarkUp)).getParamDbl("kH");
		kineticParam[paramIndex+2] = (new XMLParser(defMarkUp)).getParamDbl("nH");
		kineticParam[paramIndex+3] = (new XMLParser(defMarkUp)).getParamDbl("Ea");
		kineticParam[paramIndex+4] = (new XMLParser(defMarkUp)).getParamDbl("acidMolarMass");
		kineticParam[paramIndex+51] = (new XMLParser(defMarkUp)).getParamDbl("ironMolarMass");
		nParam = 6;
	}

	@Override
	public double kineticValue(double solute)
	{
		Double r = rate(solute, _kH, _nH, _Ea, _acidMolarMass, _ironMolarMass, _T);
		return r;
	}

	@Override
	public double kineticDiff(double solute)
	{
		Double dR = diffRate(solute, _kH, _nH, _Ea, _acidMolarMass, _ironMolarMass, _T);
		return dR;
	}

	@Override
	public double kineticValue(double solute, double[] paramTable, int index)
	{
		Double T  = paramTable[index];
		Double kH = paramTable[index+1]; 
		Double nH = paramTable[index+2];
		Double Ea = paramTable[index+3];
		Double acidMolarMass = paramTable[index+4];
		Double ironMolarMass = paramTable[index+5];

		Double r = rate(solute, kH, nH, Ea, acidMolarMass, ironMolarMass, T);

		return r;
	}

	@Override
	public double kineticDiff(double solute, double[] paramTable, int index)
	{
		Double T  = paramTable[index];
		Double kH = paramTable[index+1]; 
		Double nH = paramTable[index+2];
		Double Ea = paramTable[index+3];
		Double acidMolarMass = paramTable[index+4];
		Double ironMolarMass = paramTable[index+5];

		Double dR = diffRate(solute, kH, nH, Ea, acidMolarMass, ironMolarMass, T);

		return dR;
	}

	private Double rate(Double solute, Double kH, Double nH, Double Ea,
						Double acidMolarMass, Double ironMolarMass, Double T)
	{
		Double rate;
		if(solute > 0)
		{
			Double hPlus = hPlusConc(solute, acidMolarMass);
			rate = kH * Math.pow(hPlus, nH)
					* Math.exp(-Ea/(_R * T))
					* ironMolarMass;
 		}
		else
		{
			rate = 0.0;
		}
		return rate;
	}

	private Double diffRate(Double solute, Double kH, Double nH, Double Ea,
							Double acidMolarMass, Double ironMolarMass, Double T)
	{
		Double dR;
		if(solute > 0)
		{
			Double hPlus = hPlusConc(solute, acidMolarMass);
			dR = kH * nH * Math.pow(hPlus, nH-1)
					* Math.exp(-Ea/(_R * T))
					* ironMolarMass;
		}
		else
		{
			dR = 0.0;
		}
		return dR;
	}

	private Double hPlusConc(Double solute, Double acidMolarMass)
	{
		return solute / (acidMolarMass * _standardConc) + 10e-7;
	}
}

