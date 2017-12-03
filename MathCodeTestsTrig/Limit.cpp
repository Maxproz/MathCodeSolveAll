#include "Limit.h"

//
//double EvaluateFuncLimit(const RationalFunction<LinearFunction, LinearFunction>& InFunction)
//{
//
//
//
//}
//
//inline double EvaluateFuncLimit(const RationalFunction<QuadraticFunction, LinearFunction>& InFunction)
//{
//
//}
//
//inline double EvaluateFuncLimit(const RationalFunction<ConstantFunction, QuadraticFunction>& InFunction)
//{
//
//
//}
//
//inline double EvaluateFuncLimit(QuadraticFunction& InFunction)
//{
//
//
//}
//
//inline double EvaluateFuncLimit(const LinearFunction& InFunction)
//{
//
//	
//}
//
//inline double EvaluateFuncLimit(const CubicFunction& InFunction)
//{
//
//	
//}
//
//// TODO IMPORTANT: How do I handle the differences in power variables here?
////template <int Power>
//// TODO: Important I should move the EvaluateFuncLimit functions to the cpp file
////
//inline double EvaluateFuncLimit(const TrigometricFunction<MPSIN>& InFunction)
//{
//
//	
//}
//
////template <int Power>
//inline double EvaluateFuncLimit(const TrigometricFunction<MPCOS>& InFunction)
//{
//
//
//}
//
////template <int Power>
//inline double EvaluateFuncLimit(const TrigometricFunction<MPTAN>& InFunction)
//{
//
//
//}


//template <typename Function>
//double Limit<Function>::EvaluateRootFuncLimit(const RootFunction & InRootFunc)
//{
//	bool bDoesNegDirLimitExist = true;
//
//
//	double Limit = InRootFunc(m_a);
//
//	double DomainStart = InRootFunc.GetStartingDomainNum();
//
//
//
//	// TODO: Find a way to use the starting domain number to trigger a bool that will display 
//	// Limit Does not Exist for the side that DNE
//	// UPDATE: sort of did it needs more polishing though
//
//	std::vector<std::pair<double, double>> PosDirVec;
//	// pos direction
//	std::pair<double, double> PosPairFirst;
//	PosPairFirst.first = m_a + 0.1;
//	PosPairFirst.second = InRootFunc(m_a + 0.1);
//
//	std::pair<double, double> PosPairSecond;
//	PosPairSecond.first = m_a + 0.01;
//	PosPairSecond.second = InRootFunc(m_a + 0.01);
//
//	std::pair<double, double> PosPairThird;
//	PosPairThird.first = m_a + 0.001;
//	PosPairThird.second = InRootFunc(m_a + 0.001);
//
//	std::pair<double, double> PosPairFourth;
//	PosPairFourth.first = m_a + 0.0001;
//	PosPairFourth.second = InRootFunc(m_a + 0.0001);
//
//	PosDirVec.push_back(PosPairFirst);
//	PosDirVec.push_back(PosPairSecond);
//	PosDirVec.push_back(PosPairThird);
//	PosDirVec.push_back(PosPairFourth);
//
//	std::vector<std::pair<double, double>> NegDirVec;
//	// neg direction
//	std::pair<double, double> NegPairFirst;
//	NegPairFirst.first = m_a - 0.1;
//	NegPairFirst.second = InRootFunc(m_a - 0.1);
//
//	if (std::isnan(NegPairFirst.second))
//	{
//		bDoesNegDirLimitExist = false;
//	}
//
//	std::pair<double, double> NegPairSecond;
//	NegPairSecond.first = m_a - 0.01;
//	NegPairSecond.second = InRootFunc(m_a - 0.01);
//
//	std::pair<double, double> NegPairThird;
//	NegPairThird.first = m_a - 0.001;
//	NegPairThird.second = InRootFunc(m_a - 0.001);
//
//	std::pair<double, double> NegPairFourth;
//	NegPairFourth.first = m_a - 0.0001;
//	NegPairFourth.second = InRootFunc(m_a - 0.0001);
//
//	NegDirVec.push_back(NegPairFirst);
//	NegDirVec.push_back(NegPairSecond);
//	NegDirVec.push_back(NegPairThird);
//	NegDirVec.push_back(NegPairFourth);
//
//	std::cout << "Evaluating Limit: Please Wait...\n";
//
//	for (auto & num : PosDirVec)
//	{
//
//		std::cout << std::setprecision(std::numeric_limits<double>::digits10 + 1)
//			<< std::setw(7) << num.first << " " << num.second << std::endl;
//	}
//
//	std::cout << std::endl;
//
//	for (auto & num : NegDirVec)
//	{
//		std::cout << std::setprecision(std::numeric_limits<double>::digits10 + 1)
//			<< std::setw(7) << num.first << " " << num.second << std::endl;
//	}
//
//	double TopRes = PosDirVec[3].second;
//	double BottomRes = NegDirVec[3].second;
//	/*std::cout << TopRes << std::endl;
//	std::cout << BottomRes << std::endl;*/
//
//
//	std::string TopResStr = std::to_string(TopRes);
//	// TODO: remove debug code later
//	std::cout << TopResStr << std::endl;
//
//	std::string BottomResStr = std::to_string(BottomRes);
//	// TODO: remove debug code later
//	std::cout << BottomResStr << std::endl;
//
//	double LocalPosRes = std::stod(TopResStr) * 100;
//
//	double LocalNegRes = std::stod(BottomResStr) * 100;
//
//	//std::cout << std::setprecision(std::numeric_limits<double>::digits10 + 1)
//	//	<< std::setw(7) << std::floor(LocalPosRes) << std::endl;
//
//
//
//	std::cout << std::endl;
//
//	// Are these limits infinite?
//	bool bIsPosDirPosInfinity = false;
//	bool bIsPosDirNegInfinity = false;
//	bool bIsNegDirNegInfinity = false;
//	bool bIsNegDirPosInfinity = false;
//
//
//
//	if (PosDirVec[1].second > PosDirVec[0].second)
//	{
//		if (PosDirVec[2].second > PosDirVec[1].second)
//		{
//			bIsPosDirPosInfinity = true;
//		}
//	}
//
//	if (NegDirVec[1].second < NegDirVec[0].second)
//	{
//		if (NegDirVec[2].second < NegDirVec[1].second)
//		{
//			bIsNegDirNegInfinity = true;
//		}
//	}
//
//	if (PosDirVec[1].second < PosDirVec[0].second)
//	{
//		if (PosDirVec[2].second < PosDirVec[1].second)
//		{
//			bIsPosDirNegInfinity = true;
//		}
//	}
//	if (NegDirVec[1].second > NegDirVec[0].second)
//	{
//		if (NegDirVec[2].second > NegDirVec[1].second)
//		{
//			bIsNegDirPosInfinity = true;
//		}
//	}
//
//	if (std::ceil(LocalPosRes) == std::floor(LocalNegRes))
//	{
//		std::cout << "We have a working limit result! 1\n";
//
//
//		if (bIsPosDirPosInfinity)
//		{
//			if (bIsNegDirPosInfinity)
//			{
//				std::cout << "As x approaches " << m_a << " Limit of f(x) = ";
//				std::cout << INFINITY << std::endl;
//				return INFINITY;
//			}
//		}
//
//		if (bIsPosDirNegInfinity)
//		{
//			if (bIsNegDirNegInfinity)
//			{
//				std::cout << "As x approaches " << m_a << " Limit of f(x) = ";
//				std::cout << NEGINFINITY << std::endl;
//				return NEGINFINITY;
//			}
//		}
//
//		//std::cout << std::setprecision(std::numeric_limits<double>::digits10 + 1)
//		//	<< std::setw(7) << PosDirVec[3].first << " " << PosDirVec[3].second << std::endl;
//
//		//std::cout << "Limit: " << std::ceil(LocalPosRes) / 100 << "\n";
//
//		//std::cout << std::setprecision(std::numeric_limits<double>::digits10 + 1)
//		//	<< std::setw(7) << NegDirVec[3].first << " " << NegDirVec[3].second << std::endl;
//
//
//		//	// return either
//		//return std::ceil(LocalPosRes) / 100;
//	}
//
//	// if you have two results that are returning negatives
//	// you need to do some sort of swap with floor / ceil
//	if (std::floor(LocalPosRes) == std::ceil(LocalNegRes))
//	{
//		std::cout << "We have a working limit result! 2\n";
//
//		std::cout << std::setprecision(std::numeric_limits<double>::digits10 + 1)
//			<< std::setw(7) << PosDirVec[3].first << " " << PosDirVec[3].second << std::endl;
//
//		std::cout << "Limit: " << std::floor(LocalPosRes) / 100 << "\n";
//
//		std::cout << std::setprecision(std::numeric_limits<double>::digits10 + 1)
//			<< std::setw(7) << NegDirVec[3].first << " " << NegDirVec[3].second << std::endl;
//
//
//		// return either
//		return std::floor(LocalPosRes) / 100;
//	}
//
//
//	std::cout << "Limit: DNE (does not exist): \n\n";
//
//	std::cout << "As x approaches " << m_a << " from the positive direction f(x) = ";
//
//	if (bIsPosDirPosInfinity)
//	{
//		std::cout << INFINITY << std::endl;
//	}
//	else
//	{
//		std::cout << std::setprecision(std::numeric_limits<double>::digits10 + 1)
//			<< std::setw(7) << PosDirVec[3].first << " " << PosDirVec[3].second << std::endl;
//	}
//
//	if (bDoesNegDirLimitExist == false)
//	{
//		std::cout << "As x approaches " << m_a << " from the negative direction f(x) = DOES NOT EXIST";
//	}
//	else
//	{
//		// The limit Does Exist
//		std::cout << "As x approaches " << m_a << " from the negative direction f(x) = ";
//
//		if (bIsNegDirNegInfinity)
//		{
//			std::cout << NEGINFINITY << std::endl;
//		}
//		else
//		{
//			std::cout << std::setprecision(std::numeric_limits<double>::digits10 + 1)
//				<< std::setw(7) << NegDirVec[3].first << " " << NegDirVec[3].second << std::endl;
//		}
//	}
//
//
//
//	std::cout << std::endl;
//
//	// return 0.0???
//	return 0.0;
//}
//
