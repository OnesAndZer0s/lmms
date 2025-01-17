/*
 * ControllerDialog.cpp - per-controller-specific view for changing a
 *                        controller's settings
 *
 * Copyright (c) 2008 Paul Giblock <drfaygo/at/gmail.com>
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

#include <QCloseEvent>

#include "ControllerDialog.h"
#include "Controller.h"
#include "GuiApplication.h"
#include "MainWindow.h"


ControllerDialog::ControllerDialog( Controller * _controller,
							QWidget * _parent ) :
	QWidget( _parent ),
	ModelView( _controller, this )
{
}



ControllerDialog::~ControllerDialog()
{
}



void ControllerDialog::closeEvent( QCloseEvent * _ce )
{
	if (windowFlags().testFlag(Qt::Window))
	{
		_ce->accept();
	}
	else if (gui->mainWindow()->workspace())
	{
		parentWidget()->hide();
		_ce->ignore();
	}
	else
	{
		hide();
		_ce->ignore();
	}
	emit closed();
}




