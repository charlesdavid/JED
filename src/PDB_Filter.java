package jed;

/**
 * JED class PDB_Filter: A File Name Filter that is used when reading PDB files.
 * The extension is set to read only *.pdb files in a given directory.
 * 
 * Copyright (C) 2012 Dr. Charles David
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 * 
 * @author Dr. Charles David
 * 
 */

import java.io.File;
import java.io.FilenameFilter;

class PDB_Filter implements FilenameFilter
{

	String filter_string = ".pdb";

	/**
	 * @return The current filter string used for file filtering.
	 */
	public String get_Filter_string()
		{
			return filter_string;
		}

	/**
	 * Sets the filter string to use for file filtering.
	 * 
	 * @param text
	 *            The filter string
	 */
	public void set_Filter_string(String text)
		{
			this.filter_string = text;
		}

	/**
	 * Returns all filenames that end with the specified filter string.
	 * 
	 * @see java.io.FilenameFilter#accept(java.io.File, java.lang.String)
	 */
	@Override
	public boolean accept(File dir, String name)
		{
			return (name.endsWith(filter_string));
		}
}
