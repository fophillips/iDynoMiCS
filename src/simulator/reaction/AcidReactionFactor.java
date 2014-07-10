package simulator.reaction;

import simulator.Simulator;
import utils.XMLParser;

/**
 * Created with IntelliJ IDEA.
 * Date: 08/07/2014
 * Time: 12:08
 *
 * @author Fred Phillips
 */
public class AcidReactionFactor extends ReactionFactor
{
	/**
	 * Surface area of unit cell
	 */
	private Double _surfaceArea;
	@Override
	public void init(Simulator aSim, XMLParser xmlRoot)
	{
		super.init(aSim, xmlRoot);
		_surfaceArea = Math.pow(aSim.soluteList[0].getResolution(), 2) * 1e-12;
	}
	@Override
	public void computeUptakeRate(double[] s, double mass, double t)
	{
		super.computeSpecificGrowthRate(s);
		for (int iSolute : _mySoluteIndex) {
			_uptakeRate[iSolute] = _surfaceArea * _specRate*_soluteYield[iSolute];
		}

		int iSolute;
		for (int i = 0; i<_soluteFactor.length; i++) {
			iSolute = _soluteFactor[i];
			if(iSolute!=-1)
				_diffUptakeRate[iSolute] = _surfaceArea * marginalDiffMu[i]*_soluteYield[iSolute];
		}
	}
}
