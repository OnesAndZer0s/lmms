/*
 * AutomationPattern.h - declaration of class AutomationPattern, which contains
 *                       all information about an automation pattern
 *
 * Copyright (c) 2008-2014 Tobias Doerffel <tobydox/at/users.sourceforge.net>
 * Copyright (c) 2006-2008 Javier Serrano Polo <jasp00/at/users.sourceforge.net>
 *
 * This file is part of LMMS - https://lmms.io
 *
 * This program is free software; you can redistribute it and/or
 * modify it under the terms of the GNU General Public
 * License as published by the Free Software Foundation; either
 * version 2 of the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public
 * License along with this program (see COPYING); if not, write to the
 * Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
 * Boston, MA 02110-1301 USA.
 *
 */

#ifndef AUTOMATION_PATTERN_H
#define AUTOMATION_PATTERN_H

#include <QMap>
#include <QPair>
#include <QPointer>

#include "AutomationNode.h"
#include "TrackContentObject.h"


class AutomationTrack;
class TimePos;



class LMMS_EXPORT AutomationPattern : public TrackContentObject
{
	Q_OBJECT
public:
	enum ProgressionTypes
	{
		DiscreteProgression,
		LinearProgression,
		CubicHermiteProgression,
		BezierProgression
	} ;

<<<<<<< HEAD
	typedef QMap<int, AutomationNode> timeMap;
	typedef QVector<QPointer<AutomatableModel>> objectVector;
=======
	typedef QMap<int, float> timeMap;
	typedef QMap<int, QPair<int, float> > controlPointTimeMap;
	typedef QVector<QPointer<AutomatableModel> > objectVector;
>>>>>>> feature/bezier

	AutomationPattern( AutomationTrack * _auto_track );
	AutomationPattern( const AutomationPattern & _pat_to_copy );
	virtual ~AutomationPattern() = default;

	bool addObject( AutomatableModel * _obj, bool _search_dup = true );

	const AutomatableModel * firstObject() const;
	const objectVector& objects() const;

	// progression-type stuff
	inline ProgressionTypes progressionType() const
	{
		return m_progressionType;
	}
	void setProgressionType( ProgressionTypes _new_progression_type );

	inline float getTension() const
	{
		return m_tension;
	}
	void setTension( QString _new_tension );

	TimePos timeMapLength() const;
	void updateLength();

	TimePos putValue(
		const TimePos & time,
		const float value,
		const bool quantPos = true,
		const bool ignoreSurroundingPoints = true
	);

	TimePos putValues(
		const TimePos & time,
		const float inValue,
		const float outValue,
		const bool quantPos = true,
		const bool ignoreSurroundingPoints = true
	);

<<<<<<< HEAD
	void removeNode(const TimePos & time);
	void removeNodes(const int tick0, const int tick1);

	void resetNodes(const int tick0, const int tick1);
=======
	TimePos putControlPoint( timeMap::const_iterator it,
						const float _value);

	TimePos putControlPoint(timeMap::const_iterator it,
						const int time, const float _value);

	TimePos putControlPoint(timeMap::const_iterator it,
					const int time, const float _value, const bool flip);

	void removeValue( const TimePos & time );
>>>>>>> feature/bezier

	void recordValue(TimePos time, float value);

	TimePos setDragValue( const TimePos & time,
				const float value,
				const bool quantPos = true,
				const bool controlKey = false );

	TimePos setControlPointDragValue( const TimePos & _time, const float _value, const int _x,
						   const bool _quant_pos = true );

	void applyDragValue();

	void flipControlPoint(bool flip);


	bool isDragging() const
	{
		return m_dragging;
	}

	inline const timeMap & getTimeMap() const
	{
		return m_timeMap;
	}

	inline timeMap & getTimeMap()
	{
		return m_timeMap;
	}

<<<<<<< HEAD
=======
	inline const timeMap & getTangents() const
	{
		return m_tangents;
	}

	inline timeMap & getTangents()
	{
		return m_tangents;
	}

	inline const controlPointTimeMap & getControlPoints() const
	{
		return m_controlPoints;
	}

	inline controlPointTimeMap & getControlPoints()
	{
		return m_controlPoints;
	}

	// This is for getting the node of the control point that is being dragged
	inline const timeMap::ConstIterator & getControlPointNode() const
	{
		return m_oldControlPointNode;
	}

	inline  timeMap::ConstIterator & getControlPointNode()
	{
		return m_oldControlPointNode;
	}

>>>>>>> feature/bezier
	inline float getMin() const
	{
		return firstObject()->minValue<float>();
	}

	inline float getMax() const
	{
		return firstObject()->maxValue<float>();
	}

	inline bool hasAutomation() const
	{
		return m_timeMap.isEmpty() == false;
	}

	float valueAt( const TimePos & _time ) const;
	float *valuesAfter( const TimePos & _time ) const;

	const QString name() const;

	// settings-management
	void saveSettings( QDomDocument & _doc, QDomElement & _parent ) override;
	void loadSettings( const QDomElement & _this ) override;

	static const QString classNodeName() { return "automationpattern"; }
	QString nodeName() const override { return classNodeName(); }

	TrackContentObjectView * createView( TrackView * _tv ) override;


	static bool isAutomated( const AutomatableModel * _m );
	static QVector<AutomationPattern *> patternsForModel( const AutomatableModel * _m );
	static AutomationPattern * globalAutomationPattern( AutomatableModel * _m );
	static void resolveAllIDs();

	bool isRecording() const { return m_isRecording; }
	void setRecording( const bool b ) { m_isRecording = b; }

	static int quantization() { return s_quantization; }
	static void setQuantization(int q) { s_quantization = q; }

	void clampControlPoints(bool clampVertical=true);

public slots:
	void clear();
	void objectDestroyed( jo_id_t );
	void flipY( int min, int max );
	void flipY();
	void flipX( int length = -1 );

private:
	void cleanObjects();
	void cleanControlPoints();
	void generateTangents();
	void generateTangents(timeMap::iterator it, int numToGenerate);
	float valueAt( timeMap::const_iterator v, int offset ) const;

	// Mutex to make methods involving automation patterns thread safe
	// Mutable so we can lock it from const objects
	mutable QMutex m_patternMutex;

	AutomationTrack * m_autoTrack;
	QVector<jo_id_t> m_idsToResolve;
	objectVector m_objects;
	timeMap m_timeMap;	// actual values
	timeMap m_oldTimeMap;	// old values for storing the values before setDragValue() is called.
<<<<<<< HEAD
=======
	timeMap m_tangents;	// slope at each point for calculating spline
	controlPointTimeMap m_controlPoints;	// control points for calculating the bezier curve
	controlPointTimeMap m_oldControlPoints;	// old values for storing the values before setDragValue() is called.
	// m_oldControlPoints is similar to m_oldTimeMap, since the control points need to be dragged as well or something
	timeMap::const_iterator m_oldControlPointNode;	// Which automation point was the control point connected to?
	bool m_controlFlip; // If the lefthand control point is grabbed, the value must be flipped around the automation point

	float m_controlPointDragOffset[2];

>>>>>>> feature/bezier
	float m_tension;
	bool m_hasAutomation;
	ProgressionTypes m_progressionType;

	bool m_dragging;
<<<<<<< HEAD
	bool m_dragKeepOutValue; // Should we keep the current dragged node's outValue?
	float m_dragOutValue; // The outValue of the dragged node's
=======
>>>>>>> feature/bezier

	bool m_isRecording;
	float m_lastRecordedValue;

	static int s_quantization;

	static const float DEFAULT_MIN_VALUE;
	static const float DEFAULT_MAX_VALUE;

	friend class AutomationPatternView;
	friend class AutomationNode;

} ;


#endif
