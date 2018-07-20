//
// Copyright (C) 2013-2018 University of Amsterdam
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//


#include "factorsmodel.h"

#include <QDebug>


QVariant FactorsModel::data(const QModelIndex &index, int role) const
{
	if(index.row() >= rowCount()) return QVariant();

	if (_labels.length() <= 0) {
		return QVariant();
	}

	qDebug() << ">>>> ---------------------";
	qDebug() << QString::number(index.row());
	qDebug() << ">>>> ---------------------";

	// Column & col=_dataSet->column(index.row());

        // return QString::fromStdString(col.name());



	if(role == NameRole)
        return _labels.at(index.row());
	else if(role == TypeRole)
		return "column";

	return QVariant();
}

QHash<int, QByteArray> FactorsModel::roleNames() const {
	static const auto roles = QHash<int, QByteArray>{
		{ NameRole,					"columnName"},
		{ TypeRole,					"type"}
	};

	return roles;
}

void FactorsModel::setDataSet(DataSet *dataSet)
{
	beginResetModel();
	_dataSet = dataSet;
	endResetModel();
}

void FactorsModel::setFactors(QStringList labels)
{

	qDebug() << "___________________________";
	qDebug() << labels;
	qDebug() << "___________________________";

	beginResetModel();
	_labels = labels;
	endResetModel();
}

void FactorsModel::refreshColumn(Column * column)
{
	// int rowChanged = _dataSet->getColumnIndex(column->name());
	int rowChanged = 0;
	emit dataChanged(index(rowChanged, 0), index(rowChanged, columnCount()));
}

void FactorsModel::datasetHeaderDataChanged(Qt::Orientation orientation, int first, int last)
{
	emit dataChanged(index(first, 0), index(last, columnCount()-1));
}
