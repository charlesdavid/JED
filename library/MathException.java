/*
 * Copyright (c) 2014, Massachusetts Institute of Technology
 * Released under the BSD 2-Clause License
 * http://opensource.org/licenses/BSD-2-Clause
 */

package bits.kde;

/**
 * @author Philip DeCamp
 */
public class MathException extends Exception {

    public MathException() {
        super();
    }

    public MathException( String message ) {
        super( message );
    }

    public MathException( Throwable t ) {
        super( t );
    }

    public MathException( String message, Throwable t ) {
        super( message, t );
    }

}
