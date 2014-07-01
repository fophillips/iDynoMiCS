package simulator.reaction.kinetic;

import org.jdom.Element;
import utils.XMLParser;

/**
 * Created with IntelliJ IDEA.
 * Date: 01/07/2014
 * Time: 11:42
 *
 *
 * (52) Ea = 47500 ± 2690(J mol−1 )
 * kˆH+ =588±558(molm−2 s−1)
 * nH+ =1.16±0.107
 * kˆOH− =0.0822±0.079(molm−2 s−1)
 * nOH− =0.16±0.0309
 * @author Fred Phillips
 */
public class AcidKinetic extends IsKineticFactor
{
	/**
	 * Terature
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
	 * Gas constant
	 */
	private Double _R = 8.314; // J K−1 mol−1
	@Override
	public void init(Element defMarkUp)
	{
		_T = (new XMLParser(defMarkUp)).getParamDbl("T");
		_kH =  (new XMLParser(defMarkUp)).getParamDbl("kH");
		_nH = (new XMLParser(defMarkUp)).getParamDbl("nH");
		_kOH = (new XMLParser(defMarkUp)).getParamDbl("kOH");
		_nOH = (new XMLParser(defMarkUp)).getParamDbl("nOH");
		_Ea = (new XMLParser(defMarkUp)).getParamDbl("Ea");
		nParam = 6;
	}

	@Override
	public void initFromAgent(Element defMarkUp, double[] kineticParam, int paramIndex)
	{
		kineticParam[paramIndex] = (new XMLParser(defMarkUp)).getParamDbl("T");
		kineticParam[paramIndex+1] =  (new XMLParser(defMarkUp)).getParamDbl("kH");
		kineticParam[paramIndex+2] = (new XMLParser(defMarkUp)).getParamDbl("nH");
		kineticParam[paramIndex+3] = (new XMLParser(defMarkUp)).getParamDbl("kOH");
		kineticParam[paramIndex+4] = (new XMLParser(defMarkUp)).getParamDbl("nOH");
		kineticParam[paramIndex+5] = (new XMLParser(defMarkUp)).getParamDbl("Ea");
	}

	@Override
	public double kineticValue(double solute)
	{
		return (_kH * Math.pow(solute, -_nH) + _kOH * Math.pow(solute, _nOH) * Math.pow(10, -14 * _nOH))
				* Math.exp(-_Ea/(_R * _T));
	}

	@Override
	public double kineticDiff(double solute)
	{
		return (-_nH * _kH * Math.pow(solute, -_nH - 1)
				+ _nOH * _kOH * Math.pow(solute, _nOH-1) * Math.pow(10, -14 * _nOH))
				* Math.exp(-_Ea/(_R * _T));
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
		
		return (kH * Math.pow(solute, -nH) + kOH * Math.pow(solute, nOH) * Math.pow(10, -14 * nOH))
				* Math.exp(-Ea/(_R * T));
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
		
		return (-nH * kH * Math.pow(solute, -nH - 1)
				+ nOH * kOH * Math.pow(solute, nOH-1) * Math.pow(10, -14 * nOH))
				* Math.exp(-Ea/(_R * T));
	}
}

