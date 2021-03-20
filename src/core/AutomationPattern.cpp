/*
 * AutomationPattern.cpp - implementation of class AutomationPattern which
 *                         holds dynamic values
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

#include "AutomationPattern.h"

#include "AutomationNode.h"
#include "AutomationPatternView.h"
#include "AutomationTrack.h"
#include "BBTrackContainer.h"
#include "LocaleHelper.h"
#include "Note.h"
#include "ProjectJournal.h"
#include "Song.h"

#include <cmath>

int AutomationPattern::s_quantization = 1;
const float AutomationPattern::DEFAULT_MIN_VALUE = 0;
const float AutomationPattern::DEFAULT_MAX_VALUE = 1;


AutomationPattern::AutomationPattern( AutomationTrack * _auto_track ) :
	TrackContentObject( _auto_track ),
	m_patternMutex(QMutex::Recursive),
	m_autoTrack( _auto_track ),
	m_objects(),
	m_tension( 1.0 ),
	m_progressionType( DiscreteProgression ),
	m_dragging( false ),
	m_isRecording( false ),
	m_lastRecordedValue( 0 )
{
	changeLength( TimePos( 1, 0 ) );
	if( getTrack() )
	{
		switch( getTrack()->trackContainer()->type() )
		{
			case TrackContainer::BBContainer:
				setAutoResize( true );
				break;

			case TrackContainer::SongContainer:
				// move down
			default:
				setAutoResize( false );
				break;
		}
	}
}




AutomationPattern::AutomationPattern( const AutomationPattern & _pat_to_copy ) :
	TrackContentObject( _pat_to_copy.m_autoTrack ),
	m_patternMutex(QMutex::Recursive),
	m_autoTrack( _pat_to_copy.m_autoTrack ),
	m_objects( _pat_to_copy.m_objects ),
	m_tension( _pat_to_copy.m_tension ),
	m_progressionType( _pat_to_copy.m_progressionType )
{
	// Locks the mutex of the copied AutomationPattern to make sure it
	// doesn't change while it's being copied
	QMutexLocker m(&_pat_to_copy.m_patternMutex);

	for( timeMap::const_iterator it = _pat_to_copy.m_timeMap.begin();
				it != _pat_to_copy.m_timeMap.end(); ++it )
	{
		// Copies the automation node (in/out values and in/out tangents)
		m_timeMap[POS(it)] = it.value();
	}
	if (!getTrack()){ return; }
	switch( getTrack()->trackContainer()->type() )
	{
		case TrackContainer::BBContainer:
			setAutoResize( true );
			break;

		case TrackContainer::SongContainer:
			// move down
		default:
			setAutoResize( false );
			break;
	}
}

bool AutomationPattern::addObject( AutomatableModel * _obj, bool _search_dup )
{
	QMutexLocker m(&m_patternMutex);

	if( _search_dup && m_objects.contains(_obj) )
	{
		return false;
	}

	// the automation track is unconnected and there is nothing in the track
	if( m_objects.isEmpty() && hasAutomation() == false )
	{
		// then initialize first value
		putValue( TimePos(0), _obj->inverseScaledValue( _obj->value<float>() ), false );
	}

	m_objects += _obj;

	connect( _obj, SIGNAL( destroyed( jo_id_t ) ),
			this, SLOT( objectDestroyed( jo_id_t ) ),
						Qt::DirectConnection );

	emit dataChanged();

	return true;
}




void AutomationPattern::setProgressionType(
					ProgressionTypes _new_progression_type )
{
	QMutexLocker m(&m_patternMutex);

	if ( _new_progression_type == DiscreteProgression ||
		_new_progression_type == LinearProgression ||
		_new_progression_type == CubicHermiteProgression ||
		_new_progression_type == BezierProgression )
	{
		m_progressionType = _new_progression_type;
		emit dataChanged();
	}
}




void AutomationPattern::setTension( QString _new_tension )
{
	QMutexLocker m(&m_patternMutex);

	bool ok;
	float nt = LocaleHelper::toFloat(_new_tension, & ok);

	if( ok && nt > -0.01 && nt < 1.01 )
	{
		m_tension = nt;
	}
}




const AutomatableModel * AutomationPattern::firstObject() const
{
	QMutexLocker m(&m_patternMutex);

	AutomatableModel* model;
	if (!m_objects.isEmpty() && (model = m_objects.first()) != nullptr)
	{
		return model;
	}

	static FloatModel fm(0, DEFAULT_MIN_VALUE, DEFAULT_MAX_VALUE, 0.001);
	return &fm;
}

const AutomationPattern::objectVector& AutomationPattern::objects() const
{
	QMutexLocker m(&m_patternMutex);

	return m_objects;
}




TimePos AutomationPattern::timeMapLength() const
{
	QMutexLocker m(&m_patternMutex);

	TimePos one_bar = TimePos(1, 0);
	if (m_timeMap.isEmpty()) { return one_bar; }

	timeMap::const_iterator it = m_timeMap.end();
	tick_t last_tick = static_cast<tick_t>(POS(it - 1));
	// if last_tick is 0 (single item at tick 0)
	// return length as a whole bar to prevent disappearing TCO
	if (last_tick == 0) { return one_bar; }

	return TimePos(last_tick);
}




void AutomationPattern::updateLength()
{
	// Do not resize down in case user manually extended up
	changeLength(qMax(length(), timeMapLength()));
}




<<<<<<< HEAD
/**
 * @brief Puts an automation node on the timeMap with the given value.
 *        The inValue and outValue of the created node will be the same.
 * @param TimePos time to add the node to
 * @param Float inValue and outValue of the node
 * @param Boolean True to quantize the position (defaults to true)
 * @param Boolean True to ignore unquantized surrounding nodes (defaults to true)
 * @return TimePos of the recently added automation node
 */
TimePos AutomationPattern::putValue(
	const TimePos & time,
	const float value,
	const bool quantPos,
	const bool ignoreSurroundingPoints
)
=======
TimePos AutomationPattern::putControlPoint( timeMap::const_iterator it,
					const int time, const float _value, const bool flip )
{
	if (flip)
	{
		putControlPoint( it, 2 * it.key() - time, 2 * it.value() - _value );
	}
	else
	{
		putControlPoint( it, time, _value );
	}
	return it.key();
}




/* If we are only given the value and automation point
	then figure out where to put the control point */
TimePos AutomationPattern::putControlPoint(timeMap::const_iterator it,
							const float _value)
{
	// Insert control point at the automation point
	return putControlPoint( it, (float)it.key() + 50, _value );
}




TimePos AutomationPattern::putControlPoint(timeMap::const_iterator it,
						const int time, const float _value)
{
	m_controlPoints.remove( it.key() );
	m_controlPoints[it.key()] = {time, _value};
	clampControlPoints();
	return it.key();
}




TimePos AutomationPattern::putValue( const TimePos & time,
					const float value,
					const bool quantPos,
					const bool ignoreSurroundingPoints )
>>>>>>> feature/bezier
{
	QMutexLocker m(&m_patternMutex);

	cleanObjects();

	TimePos newTime = quantPos ? Note::quantized(time, quantization()) : time;

	// Create a node or replace the existing one on newTime
	m_timeMap[newTime] = AutomationNode(this, value, newTime);

	timeMap::iterator it = m_timeMap.find(newTime);

	// Remove control points that are covered by the new points
	// quantization value. Control Key to override
	if (!ignoreSurroundingPoints)
	{
		// We need to check that to avoid removing nodes from
		// newTime + 1 to newTime (removing the node we are adding)
		if (quantization() > 1)
		{
			// Remove nodes between the quantization points, them not
			// being included
			removeNodes(newTime + 1, newTime + quantization() - 1);
		}
	}
<<<<<<< HEAD
	if (it != m_timeMap.begin()) { --it; }
	generateTangents(it, 3);

	updateLength();

	emit dataChanged();

	return newTime;
}




/**
 * @brief Puts an automation node on the timeMap with the given inValue
 *        and outValue.
 * @param TimePos time to add the node to
 * @param Float inValue of the node
 * @param Float outValue of the node
 * @param Boolean True to quantize the position (defaults to true)
 * @param Boolean True to ignore unquantized surrounding nodes (defaults to true)
 * @return TimePos of the recently added automation node
 */
TimePos AutomationPattern::putValues(
	const TimePos & time,
	const float inValue,
	const float outValue,
	const bool quantPos,
	const bool ignoreSurroundingPoints
)
{
	QMutexLocker m(&m_patternMutex);

	cleanObjects();

	TimePos newTime = quantPos ? Note::quantized(time, quantization()) : time;

	// Create a node or replace the existing one on newTime
	m_timeMap[newTime] = AutomationNode(this, inValue, outValue, newTime);

	timeMap::iterator it = m_timeMap.find(newTime);

	// Remove control points that are covered by the new points
	// quantization value. Control Key to override
	if (!ignoreSurroundingPoints)
=======
	putControlPoint(it, value);
	clampControlPoints();

	if( it != m_timeMap.begin() )
>>>>>>> feature/bezier
	{
		// We need to check that to avoid removing nodes from
		// newTime + 1 to newTime (removing the node we are adding)
		if (quantization() > 1)
		{
			// Remove nodes between the quantization points, them not
			// being included
			removeNodes(newTime + 1, newTime + quantization() - 1);
		}
	}
	if (it != m_timeMap.begin()) { --it; }
	generateTangents(it, 3);

	updateLength();

	emit dataChanged();

	return newTime;
}




void AutomationPattern::removeNode(const TimePos & time)
{
	QMutexLocker m(&m_patternMutex);

	cleanObjects();

	m_timeMap.remove( time );
<<<<<<< HEAD
	timeMap::iterator it = m_timeMap.lowerBound(time);
=======
	m_tangents.remove( time );
	m_controlPoints.remove( time );
	timeMap::const_iterator it = m_timeMap.lowerBound( time );
>>>>>>> feature/bezier
	if( it != m_timeMap.begin() )
	{
		--it;
	}
	generateTangents(it, 3);

	updateLength();

	emit dataChanged();
}




/**
 * @brief Removes all automation nodes between the given ticks
 * @param Int first tick of the range
 * @param Int second tick of the range
 */
void AutomationPattern::removeNodes(const int tick0, const int tick1)
{
	if (tick0 == tick1)
	{
		removeNode(TimePos(tick0));
		return;
	}

	TimePos start = TimePos(qMin(tick0, tick1));
	TimePos end = TimePos(qMax(tick0, tick1));

	// Make a list of TimePos with nodes to be removed
	// because we can't simply remove the nodes from
	// the timeMap while we are iterating it.
	QVector<TimePos> nodesToRemove;

	for (auto it = m_timeMap.lowerBound(start), endIt = m_timeMap.upperBound(end); it != endIt; ++it)
	{
		nodesToRemove.append(POS(it));
	}

	for (auto node: nodesToRemove)
	{
		removeNode(node);
	}
}




/**
 * @brief Resets the outValues of all automation nodes between the given ticks
 * @param Int first tick of the range
 * @param Int second tick of the range
 */
void AutomationPattern::resetNodes(const int tick0, const int tick1)
{
	if (tick0 == tick1)
	{
		auto it = m_timeMap.find(TimePos(tick0));
		if (it != m_timeMap.end()) { it.value().resetOutValue(); }
		return;
	}

	TimePos start = TimePos(qMin(tick0, tick1));
	TimePos end = TimePos(qMax(tick0, tick1));

	for (auto it = m_timeMap.lowerBound(start), endIt = m_timeMap.upperBound(end); it != endIt; ++it)
	{
		it.value().resetOutValue();
	}
}




void AutomationPattern::recordValue(TimePos time, float value)
{
	QMutexLocker m(&m_patternMutex);

	if( value != m_lastRecordedValue )
	{
		putValue( time, value, true );
		m_lastRecordedValue = value;
	}
	else if( valueAt( time ) != value )
	{
		removeNode(time);
	}
}




/**
 * @brief Set the position of the point that is being dragged.
 *        Calling this function will also automatically set m_dragging to true.
 *        When applyDragValue() is called, m_dragging is set back to false.
 * @param TimePos of the node being dragged
 * @param Float with the value to assign to the point being dragged
 * @param Boolean. True to snip x position
 * @param Boolean. True to ignore unquantized surrounding nodes
 * @return TimePos with current time of the dragged value
 */
TimePos AutomationPattern::setDragValue(
	const TimePos & time,
	const float value,
	const bool quantPos,
	const bool controlKey
)
{
<<<<<<< HEAD
	QMutexLocker m(&m_patternMutex);

	if (m_dragging == false)
	{
		TimePos newTime = quantPos ? Note::quantized(time, quantization()) : time;

		// We will keep the same outValue only if it's different from the
		// inValue
		m_dragKeepOutValue = false;

		// Check if we already have a node on the position we are dragging
		// and if we do, store the outValue so the discrete jump can be kept
		timeMap::iterator it = m_timeMap.find(newTime);
		if (it != m_timeMap.end())
		{
			if (OFFSET(it) != 0)
			{
				m_dragKeepOutValue = true;
				m_dragOutValue = OUTVAL(it);
			}
		}

		this->removeNode(newTime);
=======
	//cleanControlPoints();
	if( m_dragging == false )
	{
		TimePos newTime = quantPos  ?
				Note::quantized( time, quantization() ) :
							time;

		if ( m_timeMap.contains( newTime ) )
		{
			// Set the offset for the control point, so it gets dragged around with the automation point
			m_controlPointDragOffset[0] = (float)m_controlPoints[newTime].first - (float)newTime;
			m_controlPointDragOffset[1] = m_controlPoints[newTime].second - m_timeMap[newTime];
		}
		else
		{
			m_controlPointDragOffset[0] = 50;
			m_controlPointDragOffset[1] = 0;
		}

		this->removeValue( newTime );
>>>>>>> feature/bezier
		m_oldTimeMap = m_timeMap;
		m_oldControlPoints = m_controlPoints;
		m_dragging = true;
	}

	//Restore to the state before it the point were being dragged
	m_timeMap = m_oldTimeMap;
	m_controlPoints = m_oldControlPoints;

	generateTangents();

	if (m_dragKeepOutValue)
	{
		return this->putValues(time, value, m_dragOutValue, quantPos, controlKey);
	}

<<<<<<< HEAD
	return this->putValue(time, value, quantPos, controlKey);
=======
	// Put the new automation point down
	TimePos returnValue = this->putValue( time, value, quantPos, controlKey );
	// Put a new control point down at an offset
	m_controlPoints.remove( returnValue );
	putControlPoint(m_timeMap.find( returnValue ), (float)returnValue + m_controlPointDragOffset[0],
			value + m_controlPointDragOffset[1]);
	clampControlPoints();

	return returnValue;
}





TimePos AutomationPattern::setControlPointDragValue( const TimePos & _time, const float _value, const int _x,
					   const bool _quant_pos)
{
	if( m_dragging == false )
	{
		TimePos newTime = _quant_pos  ?
					Note::quantized( _time, quantization() ) :
					_time;
		m_controlPoints.remove( newTime );
		m_oldControlPointNode = m_timeMap.find( newTime );
		m_dragging = true;
	}

	return this->putControlPoint(m_oldControlPointNode, _x, _value, m_controlFlip);
}




/**
 * @breif If the control point grabbed is on the left of the automation point,
 *        be flipped in order to get the control points actual location.
 * @param should the value be flipped or not
 */
void AutomationPattern::flipControlPoint(bool flip)
{
	m_controlFlip = flip;
>>>>>>> feature/bezier
}




/**
 * @brief After the point is dragged, this function is called to apply the change.
 */
void AutomationPattern::applyDragValue()
{
	QMutexLocker m(&m_patternMutex);

	m_dragging = false;
}




float AutomationPattern::valueAt( const TimePos & _time ) const
{
	QMutexLocker m(&m_patternMutex);

	if( m_timeMap.isEmpty() )
	{
		return 0;
	}

	// If we have a node at that time, just return its value
	if (m_timeMap.contains(_time))
	{
		// When the time is exactly the node's time, we want the inValue
		return m_timeMap[_time].getInValue();
	}

	// lowerBound returns next value with equal or greater key. Since we already
	// checked if the key contains a node, we know the returned node has a greater
	// key than _time. Therefore we take the previous element to calculate the current value
	timeMap::const_iterator v = m_timeMap.lowerBound(_time);

	if( v == m_timeMap.begin() )
	{
		return 0;
	}
	if( v == m_timeMap.end() )
	{
		// When the time is after the last node, we want the outValue of it
		return OUTVAL(v - 1);
	}

	return valueAt(v - 1, _time - POS(v - 1));
}




// This method will get the value at an offset from a node, so we use the outValue of
// that node and the inValue of the next node for the calculations.
float AutomationPattern::valueAt( timeMap::const_iterator v, int offset ) const
{
	QMutexLocker m(&m_patternMutex);

	// We never use it with offset 0, but doesn't hurt to return a correct
	// value if we do
	if (offset == 0) { return INVAL(v); }

	if (m_progressionType == DiscreteProgression)
	{
		return OUTVAL(v);
	}
	else if( m_progressionType == LinearProgression )
	{
		float slope =
			(INVAL(v + 1) - OUTVAL(v))
			/ (POS(v + 1) - POS(v));

		return OUTVAL(v) + offset * slope;
	}
	else if( m_progressionType == CubicHermiteProgression )
	{
		// Implements a Cubic Hermite spline as explained at:
		// http://en.wikipedia.org/wiki/Cubic_Hermite_spline#Unit_interval_.280.2C_1.29
		//
		// Note that we are not interpolating a 2 dimensional point over
		// time as the article describes.  We are interpolating a single
		// value: y.  To make this work we map the values of x that this
		// segment spans to values of t for t = 0.0 -> 1.0 and scale the
		// tangents _m1 and _m2
		int numValues = (POS(v + 1) - POS(v));
		float t = (float) offset / (float) numValues;
		float m1 = OUTTAN(v) * numValues * m_tension;
		float m2 = INTAN(v + 1) * numValues * m_tension;

		auto t2 = pow(t, 2);
		auto t3 = pow(t, 3);
		return (2 * t3 - 3 * t2 + 1) * OUTVAL(v)
			+ (t3 - 2 * t2 + t) * m1
			+ (-2 * t3 + 3 * t2) * INVAL(v + 1)
			+ (t3 - t2) * m2;
	}
	else /* BezierProgression */
	{

		/* Formula goes as such:
			Automation points: 		(x0, y0), (x3, y3)
			Relative control points: 	(x1, y1), (x2, y2)
			Where the control points are BETWEEN the automation points.
				(This isn't the case in this program, so the second control point must be "flipped"
				around its automation point)

			x = ( (1-t)^3 * x0 ) + ( 3 * (1-t)^2 * t * x1 ) + ( 3 * (1-t) * t^2 * x2 ) + ( t^3 * x3 )
			y = ( (1-t)^3 * y0 ) + ( 3 * (1-t)^2 * t * y1 ) + ( 3 * (1-t) * t^2 * y2 ) + ( t^3 * y3 )

			0 <= t <= 1
		*/

		int numValues = (v+1).key() - v.key();

		// The x values are essentially twice the distance from their control points
		// to make up for their range being limited.
		int targetX1 = ( m_controlPoints[v.key()].first - v.key() ) * 2;
		int targetX2 = ( 3 * (v+1).key() - 2 * m_controlPoints[(v+1).key()].first - v.key() );
		// The y values are the actual y values. Maybe this should be doubled,
		// but it doesn't seem necessary to me.
		float targetY1 = m_controlPoints[v.key()].second;
		float targetY2 = 2*(v+1).value() - m_controlPoints[(v+1).key()].second;

		// To find the y value on the curve at a certain x, we first have to find the t (between 0 and 1) that gives the x
		float t = 0;
		float x = 3 * pow((1-t), 2) * t * targetX1 + 3 * (1-t) * pow(t, 2) * targetX2 + pow(t, 3) * numValues;
		while (offset > x)
		{
			t += 0.05;
			x = 3 * pow((1-t), 2) * t * targetX1 + 3 * (1-t) * pow(t, 2) * targetX2 + pow(t, 3) * numValues;
		}

		float ratio = x;
		float y1 = pow((1-t),3) * v.value() + 3 * pow((1-t),2) * t * targetY1 +
				3 * (1-t) * pow(t,2) * targetY2 + pow(t,3) * (v+1).value();
		t -= 0.05;
		float y2 = pow((1-t),3) * v.value() + 3 * pow((1-t),2) * t * targetY1 +
				3 * (1-t) * pow(t,2) * targetY2 + pow(t,3) * (v+1).value();
		x = 3 * pow((1-t), 2) * t * targetX1 + 3 * (1-t) * pow(t, 2) * targetX2 + pow(t, 3) * numValues;

		// Ratio is how we get the linear extrapolation between points
		// We have to get the ratio of how close this point is to its left compared to right
		ratio = (offset - x) / (ratio - x);
		return (ratio)*y1 + (1-ratio)*y2;
	}
}




float *AutomationPattern::valuesAfter( const TimePos & _time ) const
{
	QMutexLocker m(&m_patternMutex);

	timeMap::const_iterator v = m_timeMap.lowerBound(_time);
	if( v == m_timeMap.end() || (v+1) == m_timeMap.end() )
	{
		return NULL;
	}

	int numValues = POS(v + 1) - POS(v);
	float *ret = new float[numValues];

	for( int i = 0; i < numValues; i++ )
	{
		ret[i] = valueAt( v, i );
	}

	return ret;
}




void AutomationPattern::flipY(int min, int max)
{
	QMutexLocker m(&m_patternMutex);

	bool changedTimeMap = false;

	for (auto it = m_timeMap.begin(); it != m_timeMap.end(); ++it)
	{
		// Get distance from IN/OUT values to max value
		float inValDist = max - INVAL(it);
		float outValDist = max - OUTVAL(it);

		// To flip, that will be the new distance between
		// the IN/OUT values and the min value
		it.value().setInValue(min + inValDist);
		it.value().setOutValue(min + outValDist);

<<<<<<< HEAD
		changedTimeMap = true;
=======
		if ( min < 0 )
		{
			tempValue = valueAt( ( iterate + i ).key() ) * -1;
			putValue( TimePos( (iterate + i).key() ) , tempValue, false);
			tempValue = m_controlPoints[(iterate + i).key()].second * -1;
			m_controlPoints[(iterate + i).key()].second = tempValue;
		}
		else
		{
			tempValue = max - valueAt( ( iterate + i ).key() );
			putValue( TimePos( (iterate + i).key() ) , tempValue, false);
			tempValue = max - m_controlPoints[(iterate + i).key()].second;
			m_controlPoints[(iterate + i).key()].second = tempValue;
		}
>>>>>>> feature/bezier
	}

	if (changedTimeMap)
	{
		generateTangents();
		emit dataChanged();
	}
}




void AutomationPattern::flipY()
{
	flipY(getMin(), getMax());
}




void AutomationPattern::flipX(int length)
{
	QMutexLocker m(&m_patternMutex);

	timeMap::const_iterator it = m_timeMap.lowerBound(0);

	if (it == m_timeMap.end()) { return; }

	// Temporary map where we will store the flipped version
	// of our pattern
	timeMap tempMap;
	controlPointTimeMap tempControlPoints;

	float tempValue = 0;
	float tempOutValue = 0;

	// We know the QMap isn't empty, making this safe:
	float realLength = m_timeMap.lastKey();

	// If we have a positive length, we want to flip the area covered by that
	// length, even if it goes beyond the pattern. A negative length means that
	// we just want to flip the nodes we have
	if (length >= 0 && length != realLength)
	{
		// If length to be flipped is bigger than the real length
		if (realLength < length)
		{
			// We are flipping an area that goes beyond the last node. So we add a node to the
			// beginning of the flipped timeMap representing the value of the end of the area
			tempValue = valueAt(length);
			tempMap[0] = AutomationNode(this, tempValue, 0);

			// Now flip the nodes we have in relation to the length
			do
			{
<<<<<<< HEAD
				// We swap the inValue and outValue when flipping horizontally
				tempValue = OUTVAL(it);
				tempOutValue = INVAL(it);
				TimePos newTime = TimePos(length - POS(it));

				tempMap[newTime] = AutomationNode(this, tempValue, tempOutValue, newTime);

				++it;
			} while (it != m_timeMap.end());
=======
				tempValue = valueAt( ( iterate + i ).key() );
				TimePos newTime = TimePos( length - ( iterate + i ).key() );

				int newControlPointX = -( iterate + i ).key() + m_controlPoints[( iterate + i ).key()].first + newTime;
				tempControlPoints[newTime] = {newControlPointX,
						2*tempValue - m_controlPoints[( iterate + i ).key()].second};

				tempMap[newTime] = tempValue;
			}
>>>>>>> feature/bezier
		}
		else // If the length to be flipped is smaller than the real length
		{
			do
			{
				TimePos newTime;

<<<<<<< HEAD
				// Only flips the length to be flipped and keep the remaining values in place
				// We also only swap the inValue and outValue if we are flipping the node
				if (POS(it) <= length)
=======
				int newControlPointX = -( iterate + i ).key() + m_controlPoints[( iterate + i ).key()].first + newTime;
				tempControlPoints[newTime] = {newControlPointX,
						2*tempValue - m_controlPoints[( iterate + i ).key()].second};

				if ( ( iterate + i ).key() <= length )
>>>>>>> feature/bezier
				{
					newTime = length - POS(it);
					tempValue = OUTVAL(it);
					tempOutValue = INVAL(it);
				}
				else
				{
					newTime = POS(it);
					tempValue = INVAL(it);
					tempOutValue = OUTVAL(it);
				}

				tempMap[newTime] = AutomationNode(this, tempValue, tempOutValue, newTime);

				++it;
			} while (it != m_timeMap.end());
		}
	}
	else // Length to be flipped is the same as the real length
	{
		do
		{
<<<<<<< HEAD
			// Swap the inValue and outValue
			tempValue = OUTVAL(it);
			tempOutValue = INVAL(it);

			TimePos newTime = TimePos(realLength - POS(it));
			tempMap[newTime] = AutomationNode(this, tempValue, tempOutValue, newTime);

			++it;
		} while (it != m_timeMap.end());
=======
			tempValue = valueAt( ( iterate + i ).key() );
			cleanObjects();
			TimePos newTime = TimePos( realLength - ( iterate + i ).key() );
			int newControlPointX = -( iterate + i ).key() + m_controlPoints[( iterate + i ).key()].first + newTime;

			tempMap[newTime] = tempValue;
			tempControlPoints[newTime] = {newControlPointX,
					2*tempValue - m_controlPoints[( iterate + i ).key()].second};
		}
>>>>>>> feature/bezier
	}

	m_timeMap.clear();
	m_controlPoints.clear();

	m_timeMap = tempMap;
	m_controlPoints = tempControlPoints;

	cleanObjects();

	generateTangents();
	emit dataChanged();
}




void AutomationPattern::saveSettings( QDomDocument & _doc, QDomElement & _this )
{
	QMutexLocker m(&m_patternMutex);

	_this.setAttribute( "pos", startPosition() );
	_this.setAttribute( "len", length() );
	_this.setAttribute( "name", name() );
	_this.setAttribute( "prog", QString::number( progressionType() ) );
	_this.setAttribute( "tens", QString::number( getTension() ) );
	_this.setAttribute( "mute", QString::number( isMuted() ) );
	
	if( usesCustomClipColor() )
	{
		_this.setAttribute( "color", color().name() );
	}

	for( timeMap::const_iterator it = m_timeMap.begin();
						it != m_timeMap.end(); ++it )
	{
		QDomElement element = _doc.createElement( "time" );
		element.setAttribute("pos", POS(it));
		element.setAttribute("value", INVAL(it));
		element.setAttribute("outValue", OUTVAL(it));
		_this.appendChild( element );
	}

	for( controlPointTimeMap::const_iterator it = m_controlPoints.begin();
						it != m_controlPoints.end(); ++it )
	{
		QDomElement element = _doc.createElement( "ctrlpnt" );
		element.setAttribute( "pos", it.key() );
		element.setAttribute( "value1", it.value().first );
		element.setAttribute( "value2", it.value().second );
		_this.appendChild( element );
	}

	for( objectVector::const_iterator it = m_objects.begin();
						it != m_objects.end(); ++it )
	{
		if( *it )
		{
			QDomElement element = _doc.createElement( "object" );
			element.setAttribute( "id",
				ProjectJournal::idToSave( ( *it )->id() ) );
			_this.appendChild( element );
		}
	}
}




void AutomationPattern::loadSettings( const QDomElement & _this )
{
	QMutexLocker m(&m_patternMutex);

	clear();

	movePosition( _this.attribute( "pos" ).toInt() );
	setName( _this.attribute( "name" ) );
	setProgressionType( static_cast<ProgressionTypes>( _this.attribute(
							"prog" ).toInt() ) );
	setTension( _this.attribute( "tens" ) );
	setMuted(_this.attribute( "mute", QString::number( false ) ).toInt() );

	for( QDomNode node = _this.firstChild(); !node.isNull();
						node = node.nextSibling() )
	{
		QDomElement element = node.toElement();
		if( element.isNull()  )
		{
			continue;
		}
		if( element.tagName() == "time" )
		{
			int timeMapPos = element.attribute("pos").toInt();
			float timeMapInValue = LocaleHelper::toFloat(element.attribute("value"));
			float timeMapOutValue = LocaleHelper::toFloat(element.attribute("outValue"));

			m_timeMap[timeMapPos] = AutomationNode(this, timeMapInValue, timeMapOutValue, timeMapPos);
		}
		else if( element.tagName() == "ctrlpnt" )
		{
			m_controlPoints[element.attribute( "pos" ).toInt()] = {element.attribute( "value1" ).toInt(),
					element.attribute( "value2" ).toFloat()};
		}
		else if( element.tagName() == "object" )
		{
			m_idsToResolve << element.attribute( "id" ).toInt();
		}
	}
	
	if( _this.hasAttribute( "color" ) )
	{
		useCustomClipColor( true );
		setColor( _this.attribute( "color" ) );
	}

	int len = _this.attribute( "len" ).toInt();
	if( len <= 0 )
	{
		// TODO: Handle with an upgrade method
		updateLength();
	}
	else
	{
		changeLength( len );
	}
	generateTangents();

	// Very important for reading older files
	cleanControlPoints();
}




const QString AutomationPattern::name() const
{
	QMutexLocker m(&m_patternMutex);

	if( !TrackContentObject::name().isEmpty() )
	{
		return TrackContentObject::name();
	}
	if( !m_objects.isEmpty() && m_objects.first() != NULL )
	{
		return m_objects.first()->fullDisplayName();
	}
	return tr( "Drag a control while pressing <%1>" ).arg(UI_CTRL_KEY);
}




void AutomationPattern::clampControlPoints(bool clampVertical)
{
	timeMap::const_iterator it;
	for (it = m_timeMap.begin(); it != m_timeMap.end(); it++)
	{
		int new_x = m_controlPoints[it.key()].first;
		float new_y = m_controlPoints[it.key()].second;
		// Clamp X positions
		// If the control point x is less than its automation point
		if ( it.key() > new_x )
		{
			new_x = it.key();
		}
		// The control point x must not pass the midpoints of its automation point and the automation points
		// its left and right
		else if ( it != m_timeMap.begin() && it.key() * 2 - new_x < ( (it-1).key() + it.key() ) / 2 )
		{
			new_x = it.key() * 2 - ( (it-1).key() + it.key() )/2;
		}
		else if ( it+1 != m_timeMap.end() && new_x > ( (it+1).key() + it.key() )/2 )
		{
			new_x = ( (it+1).key() + it.key() )/2;
		}

		if (clampVertical)
		{
			// Clamp y positions between the top and bottom of the screen
			// Clamp the right control point (keep in mind the last control point isn't clamped)
			if ( it+1 != m_timeMap.end() && new_y > getMax() )
			{
				new_y = getMax();
			}
			else if ( it+1 != m_timeMap.end() && new_y < getMin() )
			{
				new_y = getMin();
			}
			// Clamp the left control point (keep in mind the first control point isn't clamped)
			if ( it != m_timeMap.begin() && 2 * it.value() - new_y > getMax() )
			{
				new_y = 2 * it.value() - getMax();
			}
			else if ( it != m_timeMap.begin() && 2 * it.value() - new_y < getMin() )
			{
				new_y =  2 * it.value() - getMin();
			}
		}

		m_controlPoints.remove( it.key() );

		m_controlPoints[it.key()] = {new_x, new_y};
	}
}




TrackContentObjectView * AutomationPattern::createView( TrackView * _tv )
{
	QMutexLocker m(&m_patternMutex);

	return new AutomationPatternView( this, _tv );
}





bool AutomationPattern::isAutomated( const AutomatableModel * _m )
{
	TrackContainer::TrackList l;
	l += Engine::getSong()->tracks();
	l += Engine::getBBTrackContainer()->tracks();
	l += Engine::getSong()->globalAutomationTrack();

	for( TrackContainer::TrackList::ConstIterator it = l.begin(); it != l.end(); ++it )
	{
		if( ( *it )->type() == Track::AutomationTrack ||
			( *it )->type() == Track::HiddenAutomationTrack )
		{
			const Track::tcoVector & v = ( *it )->getTCOs();
			for( Track::tcoVector::ConstIterator j = v.begin(); j != v.end(); ++j )
			{
				const AutomationPattern * a = dynamic_cast<const AutomationPattern *>( *j );
				if( a && a->hasAutomation() )
				{
					for( objectVector::const_iterator k = a->m_objects.begin(); k != a->m_objects.end(); ++k )
					{
						if( *k == _m )
						{
							return true;
						}
					}
				}
			}
		}
	}
	return false;
}


/**
 * @brief returns a list of all the automation patterns that are connected to a specific model
 * @param _m the model we want to look for
 */
QVector<AutomationPattern *> AutomationPattern::patternsForModel( const AutomatableModel * _m )
{
	QVector<AutomationPattern *> patterns;
	TrackContainer::TrackList l;
	l += Engine::getSong()->tracks();
	l += Engine::getBBTrackContainer()->tracks();
	l += Engine::getSong()->globalAutomationTrack();

	// go through all tracks...
	for( TrackContainer::TrackList::ConstIterator it = l.begin(); it != l.end(); ++it )
	{
		// we want only automation tracks...
		if( ( *it )->type() == Track::AutomationTrack ||
			( *it )->type() == Track::HiddenAutomationTrack )
		{
			// get patterns in those tracks....
			const Track::tcoVector & v = ( *it )->getTCOs();
			// go through all the patterns...
			for( Track::tcoVector::ConstIterator j = v.begin(); j != v.end(); ++j )
			{
				AutomationPattern * a = dynamic_cast<AutomationPattern *>( *j );
				// check that the pattern has automation
				if( a && a->hasAutomation() )
				{
					// now check is the pattern is connected to the model we want by going through all the connections
					// of the pattern
					bool has_object = false;
					for( objectVector::const_iterator k = a->m_objects.begin(); k != a->m_objects.end(); ++k )
					{
						if( *k == _m )
						{
							has_object = true;
						}
					}
					// if the patterns is connected to the model, add it to the list
					if( has_object ) { patterns += a; }
				}
			}
		}
	}
	return patterns;
}



AutomationPattern * AutomationPattern::globalAutomationPattern(
							AutomatableModel * _m )
{
	AutomationTrack * t = Engine::getSong()->globalAutomationTrack();
	Track::tcoVector v = t->getTCOs();
	for( Track::tcoVector::const_iterator j = v.begin(); j != v.end(); ++j )
	{
		AutomationPattern * a = dynamic_cast<AutomationPattern *>( *j );
		if( a )
		{
			for( objectVector::const_iterator k = a->m_objects.begin();
												k != a->m_objects.end(); ++k )
			{
				if( *k == _m )
				{
					return a;
				}
			}
		}
	}

	AutomationPattern * a = new AutomationPattern( t );
	a->addObject( _m, false );
	return a;
}




void AutomationPattern::resolveAllIDs()
{
	TrackContainer::TrackList l = Engine::getSong()->tracks() +
				Engine::getBBTrackContainer()->tracks();
	l += Engine::getSong()->globalAutomationTrack();
	for( TrackContainer::TrackList::iterator it = l.begin();
							it != l.end(); ++it )
	{
		if( ( *it )->type() == Track::AutomationTrack ||
			 ( *it )->type() == Track::HiddenAutomationTrack )
		{
			Track::tcoVector v = ( *it )->getTCOs();
			for( Track::tcoVector::iterator j = v.begin();
							j != v.end(); ++j )
			{
				AutomationPattern * a = dynamic_cast<AutomationPattern *>( *j );
				if( a )
				{
					for( QVector<jo_id_t>::Iterator k = a->m_idsToResolve.begin();
									k != a->m_idsToResolve.end(); ++k )
					{
						JournallingObject * o = Engine::projectJournal()->
														journallingObject( *k );
						if( o && dynamic_cast<AutomatableModel *>( o ) )
						{
							a->addObject( dynamic_cast<AutomatableModel *>( o ), false );
						}
						else
						{
							// FIXME: Remove this block once the automation system gets fixed
							// This is a temporary fix for https://github.com/LMMS/lmms/issues/3781
							o = Engine::projectJournal()->journallingObject(ProjectJournal::idFromSave(*k));
							if( o && dynamic_cast<AutomatableModel *>( o ) )
							{
								a->addObject( dynamic_cast<AutomatableModel *>( o ), false );
							}
							else
							{
								// FIXME: Remove this block once the automation system gets fixed
								// This is a temporary fix for https://github.com/LMMS/lmms/issues/4781
								o = Engine::projectJournal()->journallingObject(ProjectJournal::idToSave(*k));
								if( o && dynamic_cast<AutomatableModel *>( o ) )
								{
									a->addObject( dynamic_cast<AutomatableModel *>( o ), false );
								}
							}
						}
					}
					a->m_idsToResolve.clear();
					a->dataChanged();
				}
			}
		}
	}
}




void AutomationPattern::clear()
{
	QMutexLocker m(&m_patternMutex);

	m_timeMap.clear();
<<<<<<< HEAD
=======
	m_controlPoints.clear();
	m_tangents.clear();
>>>>>>> feature/bezier

	emit dataChanged();
}




void AutomationPattern::objectDestroyed( jo_id_t _id )
{
	QMutexLocker m(&m_patternMutex);

	// TODO: distict between temporary removal (e.g. LADSPA controls
	// when switching samplerate) and real deletions because in the latter
	// case we had to remove ourselves if we're the global automation
	// pattern of the destroyed object
	m_idsToResolve += _id;

	for( objectVector::Iterator objIt = m_objects.begin();
		objIt != m_objects.end(); objIt++ )
	{
		Q_ASSERT( !(*objIt).isNull() );
		if( (*objIt)->id() == _id )
		{
			//Assign to objIt so that this loop work even break; is removed.
			objIt = m_objects.erase( objIt );
			break;
		}
	}

	emit dataChanged();
}




void AutomationPattern::cleanObjects()
{
	QMutexLocker m(&m_patternMutex);

	for( objectVector::iterator it = m_objects.begin(); it != m_objects.end(); )
	{
		if( *it )
		{
			++it;
		}
		else
		{
			it = m_objects.erase( it );
		}
	}
}




void AutomationPattern::cleanControlPoints()
{
	// If there's any control points that aren't connected to an automation point then destroy it
	for( controlPointTimeMap::iterator it = m_controlPoints.begin(); it != m_controlPoints.end(); )
	{
		if(m_timeMap.contains( (int)it.key()) )
		{
			it++;
		}
		else
		{
			it = m_controlPoints.erase( it );
		}
	}

	// If there's any automation points without a control point then insert control points
	for( timeMap::iterator it = m_timeMap.begin(); it != m_timeMap.end(); )
	{
		if(m_controlPoints.contains( (int)it.key()) )
		{
			it++;
		}
		else
		{
			m_controlPoints[it.key()] = {it.key() + 50, it.value()};
		}
	}

	clampControlPoints(false);
}




void AutomationPattern::generateTangents()
{
	generateTangents(m_timeMap.begin(), m_timeMap.size());
}




// We have two tangents, one for the left side of the node and one for the right side
// of the node (in case we have discrete value jumps in the middle of a curve).
// If the inValue and outValue of a node are the same, consequently the inTangent and
// outTangent values of the node will be the same too.
void AutomationPattern::generateTangents(timeMap::iterator it, int numToGenerate)
{
	QMutexLocker m(&m_patternMutex);

	if( m_timeMap.size() < 2 && numToGenerate > 0 )
	{
		it.value().setInTangent(0);
		it.value().setOutTangent(0);
		return;
	}

	for( int i = 0; i < numToGenerate; i++ )
	{
		if( it == m_timeMap.begin() )
		{
			// On the first node there's no curve behind it, so we will only calculate the outTangent
			// and inTangent will be set to 0.
			float tangent = (INVAL(it + 1) - OUTVAL(it)) / (POS(it + 1) - POS(it));
			it.value().setInTangent(0);
			it.value().setOutTangent(tangent);
		}
		else if( it+1 == m_timeMap.end() )
		{
			// Previously, the last value's tangent was always set to 0. That logic was kept for both tangents
			// of the last node
			it.value().setInTangent(0);
			it.value().setOutTangent(0);
			return;
		}
		else
		{
			// When we are in a node that is in the middle of two other nodes, we need to check if we
			// have a discrete jump at this node. If we do not, then we can calculate the tangents normally.
			// If we do have a discrete jump, then we have to calculate the tangents differently for each side
			// of the curve.
			// TODO: This behavior means that a very small difference between the inValue and outValue can
			// result in a big change in the curve. In the future, allowing the user to manually adjust
			// the tangents would be better.
			float inTangent;
			float outTangent;
			if (OFFSET(it) == 0)
			{
				inTangent = (INVAL(it + 1) - OUTVAL(it - 1)) / (POS(it + 1) - POS(it - 1));
				it.value().setInTangent(inTangent);
				// inTangent == outTangent in this case
				it.value().setOutTangent(inTangent);
			}
			else
			{
				// Calculate the left side of the curve
				inTangent = (INVAL(it) - OUTVAL(it - 1)) / (POS(it) - POS(it - 1));
				// Calculate the right side of the curve
				outTangent = (INVAL(it + 1) - OUTVAL(it)) / (POS(it + 1) - POS(it));
				it.value().setInTangent(inTangent);
				it.value().setOutTangent(outTangent);
			}
		}
		it++;
	}
}
