#include "labelfiltergenerator.h"

labelFilterGenerator::labelFilterGenerator(DataSetPackage *package, QObject *parent) : QObject(parent)
{
	_package = package;
}

void labelFilterGenerator::labelFilterChanged()
{

	int neededFilters = 0;

	for(Column & col : _package->dataSet->columns())
		if(labelNeedsFilter(col))
			neededFilters++;

	std::stringstream newGeneratedFilter;

	newGeneratedFilter << "genFilter <- ";

	if(neededFilters == 0)
		newGeneratedFilter << "rep(TRUE," << _package->dataSet->rowCount() << ")";
	else
	{
		bool moreThanOne = neededFilters > 1, first = true;

		if(moreThanOne)
			newGeneratedFilter << "(";


		for(Column & col : _package->dataSet->columns())
			if(labelNeedsFilter(col))
			{
				newGeneratedFilter << (first ? "" : " & \n") << generateLabelFilter(col);
				first = false;
			}

		if(moreThanOne)
			newGeneratedFilter << ")";
	}


	emit setGeneratedFilter(QString::fromStdString(newGeneratedFilter.str()));
}

bool labelFilterGenerator::labelNeedsFilter(Column & column)
{
	for(const Label & label : column.labels())
		if(!label.filterAllows())
			return true;

	return false;
}

std::string	labelFilterGenerator::generateLabelFilter(Column & column)
{
	std::string columnName = column.name();
	std::stringstream out;
	int pos = 0, neg = 0;
	bool first = true;

	for(const Label & label : column.labels())
		(label.filterAllows() ? pos : neg)++;

	bool bePositive = pos <= neg;

	out << "(";

	for(const Label & label : column.labels())
		if(label.filterAllows() == bePositive)
		{
			out << (!first ? (bePositive ? " | " : " & ") : "") << columnName << (bePositive ? " == \"" : " != \"") << label.text() << "\"";
			first = false;
		}
	out << ")";

	return out.str();
}
