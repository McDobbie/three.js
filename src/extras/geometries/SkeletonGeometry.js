/**
 * Created by mimicry on 2014/11/8.
 */

function Matrix() {}
Matrix.prototype = {

    // Returns element (i,j) of the matrix
    e: function(i,j) {
        if (i < 1 || i > this.elements.length || j < 1 || j > this.elements[0].length) { return null; }
        return this.elements[i-1][j-1];
    },

    // Returns the number of rows/columns the matrix has
    dimensions: function() {
        return {rows: this.elements.length, cols: this.elements[0].length};
    },

    // Returns the number of rows in the matrix
    rows: function() {
        return this.elements.length;
    },

    // Returns the number of columns in the matrix
    cols: function() {
        return this.elements[0].length;
    },


    // Returns a copy of the matrix
    dup: function() {
        return Matrix.create(this.elements);
    },

    // Maps the matrix to another matrix (of the same dimensions) according to the given function
    map: function(fn) {
        var els = [], ni = this.elements.length, ki = ni, i, nj, kj = this.elements[0].length, j;
        do { i = ki - ni;
            nj = kj;
            els[i] = [];
            do { j = kj - nj;
                els[i][j] = fn(this.elements[i][j], i + 1, j + 1);
            } while (--nj);
        } while (--ni);
        return Matrix.create(els);
    },

    // Returns true iff the argument has the same dimensions as the matrix
    isSameSizeAs: function(matrix) {
        var M = matrix.elements || matrix;
        if (typeof(M[0][0]) == 'undefined') { M = Matrix.create(M).elements; }
        return (this.elements.length == M.length &&
            this.elements[0].length == M[0].length);
    },

    // Returns the result of adding the argument to the matrix
    add: function(matrix) {
        var M = matrix.elements || matrix;
        if (typeof(M[0][0]) == 'undefined') { M = Matrix.create(M).elements; }
        if (!this.isSameSizeAs(M)) { return null; }
        return this.map(function(x, i, j) { return x + M[i-1][j-1]; });
    },

    // Returns the result of subtracting the argument from the matrix
    subtract: function(matrix) {
        var M = matrix.elements || matrix;
        if (typeof(M[0][0]) == 'undefined') { M = Matrix.create(M).elements; }
        if (!this.isSameSizeAs(M)) { return null; }
        return this.map(function(x, i, j) { return x - M[i-1][j-1]; });
    },

    // Returns true iff the matrix can multiply the argument from the left
    canMultiplyFromLeft: function(matrix) {
        var M = matrix.elements || matrix;
        if (typeof(M[0][0]) == 'undefined') { M = Matrix.create(M).elements; }
        // this.columns should equal matrix.rows
        return (this.elements[0].length == M.length);
    },

    // Returns the result of multiplying the matrix from the right by the argument.
    // If the argument is a scalar then just multiply all the elements. If the argument is
    // a vector, a vector is returned, which saves you having to remember calling
    // col(1) on the result.
    multiply: function(matrix) {
        if (!matrix.elements) {
            return this.map(function(x) { return x * matrix; });
        }
        var returnVector = matrix.modulus ? true : false;
        var M = matrix.elements || matrix;
        if (typeof(M[0][0]) == 'undefined') { M = Matrix.create(M).elements; }
        if (!this.canMultiplyFromLeft(M)) { return null; }
        var ni = this.elements.length, ki = ni, i, nj, kj = M[0].length, j;
        var cols = this.elements[0].length, elements = [], sum, nc, c;
        do { i = ki - ni;
            elements[i] = [];
            nj = kj;
            do { j = kj - nj;
                sum = 0;
                nc = cols;
                do { c = cols - nc;
                    sum += this.elements[i][c] * M[c][j];
                } while (--nc);
                elements[i][j] = sum;
            } while (--nj);
        } while (--ni);
        var M = Matrix.create(elements);
        return returnVector ? M.col(1) : M;
    },

    x: function(matrix) { return this.multiply(matrix); },

    // Returns a submatrix taken from the matrix
    // Argument order is: start row, start col, nrows, ncols
    // Element selection wraps if the required index is outside the matrix's bounds, so you could
    // use this to perform row/column cycling or copy-augmenting.
    minor: function(a, b, c, d) {
        var elements = [], ni = c, i, nj, j;
        var rows = this.elements.length, cols = this.elements[0].length;
        do { i = c - ni;
            elements[i] = [];
            nj = d;
            do { j = d - nj;
                elements[i][j] = this.elements[(a+i-1)%rows][(b+j-1)%cols];
            } while (--nj);
        } while (--ni);
        return Matrix.create(elements);
    },

    // Returns the transpose of the matrix
    transpose: function() {
        var rows = this.elements.length, cols = this.elements[0].length;
        var elements = [], ni = cols, i, nj, j;
        do { i = cols - ni;
            elements[i] = [];
            nj = rows;
            do { j = rows - nj;
                elements[i][j] = this.elements[j][i];
            } while (--nj);
        } while (--ni);
        return Matrix.create(elements);
    },

    // Returns true iff the matrix is square
    isSquare: function() {
        return (this.elements.length == this.elements[0].length);
    },

    // Returns the (absolute) largest element of the matrix
    max: function() {
        var m = 0, ni = this.elements.length, ki = ni, i, nj, kj = this.elements[0].length, j;
        do { i = ki - ni;
            nj = kj;
            do { j = kj - nj;
                if (Math.abs(this.elements[i][j]) > Math.abs(m)) { m = this.elements[i][j]; }
            } while (--nj);
        } while (--ni);
        return m;
    },

    // Returns the indeces of the first match found by reading row-by-row from left to right
    indexOf: function(x) {
        var index = null, ni = this.elements.length, ki = ni, i, nj, kj = this.elements[0].length, j;
        do { i = ki - ni;
            nj = kj;
            do { j = kj - nj;
                if (this.elements[i][j] == x) { return {i: i+1, j: j+1}; }
            } while (--nj);
        } while (--ni);
        return null;
    },


    // Make the matrix upper (right) triangular by Gaussian elimination.
    // This method only adds multiples of rows to other rows. No rows are
    // scaled up or switched, and the determinant is preserved.
    toRightTriangular: function() {
        var M = this.dup(), els;
        var n = this.elements.length, k = n, i, np, kp = this.elements[0].length, p;
        do { i = k - n;
            if (M.elements[i][i] == 0) {
                for (j = i + 1; j < k; j++) {
                    if (M.elements[j][i] != 0) {
                        els = []; np = kp;
                        do { p = kp - np;
                            els.push(M.elements[i][p] + M.elements[j][p]);
                        } while (--np);
                        M.elements[i] = els;
                        break;
                    }
                }
            }
            if (M.elements[i][i] != 0) {
                for (j = i + 1; j < k; j++) {
                    var multiplier = M.elements[j][i] / M.elements[i][i];
                    els = []; np = kp;
                    do { p = kp - np;
                        // Elements with column numbers up to an including the number
                        // of the row that we're subtracting can safely be set straight to
                        // zero, since that's the point of this routine and it avoids having
                        // to loop over and correct rounding errors later
                        els.push(p <= i ? 0 : M.elements[j][p] - M.elements[i][p] * multiplier);
                    } while (--np);
                    M.elements[j] = els;
                }
            }
        } while (--n);
        return M;
    },

    toUpperTriangular: function() { return this.toRightTriangular(); },

    // Returns the determinant for square matrices
    determinant: function() {
        if (!this.isSquare()) { return null; }
        var M = this.toRightTriangular();
        var det = M.elements[0][0], n = M.elements.length - 1, k = n, i;
        do { i = k - n + 1;
            det = det * M.elements[i][i];
        } while (--n);
        return det;
    },

    det: function() { return this.determinant(); },

    // Returns true iff the matrix is singular
    isSingular: function() {
        return (this.isSquare() && this.determinant() === 0);
    },

    // Returns the trace for square matrices
    trace: function() {
        if (!this.isSquare()) { return null; }
        var tr = this.elements[0][0], n = this.elements.length - 1, k = n, i;
        do { i = k - n + 1;
            tr += this.elements[i][i];
        } while (--n);
        return tr;
    },

    tr: function() { return this.trace(); },
    rk: function() { return this.rank(); },

    // Returns the result of attaching the given argument to the right-hand side of the matrix
    augment: function(matrix) {
        var M = matrix.elements || matrix;
        if (typeof(M[0][0]) == 'undefined') { M = Matrix.create(M).elements; }
        var T = this.dup(), cols = T.elements[0].length;
        var ni = T.elements.length, ki = ni, i, nj, kj = M[0].length, j;
        if (ni != M.length) { return null; }
        do { i = ki - ni;
            nj = kj;
            do { j = kj - nj;
                T.elements[i][cols + j] = M[i][j];
            } while (--nj);
        } while (--ni);
        return T;
    },
    print: function() {
        var i = this.elements.length,
            j = this.elements[0].length,
            s,
            count, count2;
        s = "[";
        for (count = 0; count < i; count += 1) {
            for (count2 = 0; count2 < j; count2 += 1) {
                if (count2 === j-1) {
                    s += " " + this.elements[count][count2] + ";";
                }else {
                    s += " " + this.elements[count][count2] + ",";
                }

            }
        }
        s += "]";

        console.log(s);
    },

    // Returns the inverse (if one exists) using Gauss-Jordan
    inverse: function() {
        if (!this.isSquare() || this.isSingular()) { return null; }
        var ni = this.elements.length, ki = ni, i, j;
        var M = this.augment(Matrix.I(ni)).toRightTriangular();
        var np, kp = M.elements[0].length, p, els, divisor;
        var inverse_elements = [], new_element;
        // Matrix is non-singular so there will be no zeros on the diagonal
        // Cycle through rows from last to first
        do { i = ni - 1;
            // First, normalise diagonal elements to 1
            els = []; np = kp;
            inverse_elements[i] = [];
            divisor = M.elements[i][i];
            do { p = kp - np;
                new_element = M.elements[i][p] / divisor;
                els.push(new_element);
                // Shuffle of the current row of the right hand side into the results
                // array as it will not be modified by later runs through this loop
                if (p >= ki) { inverse_elements[i].push(new_element); }
            } while (--np);
            M.elements[i] = els;
            // Then, subtract this row from those above it to
            // give the identity matrix on the left hand side
            for (j = 0; j < i; j++) {
                els = []; np = kp;
                do { p = kp - np;
                    els.push(M.elements[j][p] - M.elements[i][p] * M.elements[j][i]);
                } while (--np);
                M.elements[j] = els;
            }
        } while (--ni);
        return Matrix.create(inverse_elements);
    },

    inv: function() { return this.inverse(); },

    // Returns the result of rounding all the elements
    round: function() {
        return this.map(function(x) { return Math.round(x); });
    },



    // Set the matrix's elements from an array. If the argument passed
    // is a vector, the resulting matrix will be a single column.
    setElements: function(els) {
        var i, elements = els.elements || els;
        if (typeof(elements[0][0]) != 'undefined') {
            var ni = elements.length, ki = ni, nj, kj, j;
            this.elements = [];
            do { i = ki - ni;
                nj = elements[i].length; kj = nj;
                this.elements[i] = [];
                do { j = kj - nj;
                    this.elements[i][j] = elements[i][j];
                } while (--nj);
            } while(--ni);
            return this;
        }
        var n = elements.length, k = n;
        this.elements = [];
        do { i = k - n;
            this.elements.push([elements[i]]);
        } while (--n);
        return this;
    }
};

// Constructor function
Matrix.create = function(elements) {
    var M = new Matrix();
    return M.setElements(elements);
};

// Identity matrix of size n
Matrix.I = function(n) {
    var els = [], k = n, i, nj, j;
    do { i = k - n;
        els[i] = []; nj = k;
        do { j = k - nj;
            els[i][j] = (i == j) ? 1 : 0;
        } while (--nj);
    } while (--n);
    return Matrix.create(els);
};

Matrix.Zero = function(n, m) {
    var els = [], ni = n, i, nj, j;
    do { i = n - ni;
        els[i] = [];
        nj = m;
        do { j = m - nj;
            els[i][j] = 0;
        } while (--nj);
    } while (--ni);
    return Matrix.create(els);
};

THREE.SkeletonGeometry = function (geo) {
        THREE.Geometry.call(this);
        this.geo = geo || null;
        /**
         *MyEdge
         * @param i  vertex
         * @param j  vertex
         * @param x1  the third coplanar vertex of the 1st face
         * @param x2  the third coplanar vertex of the  2d face
         * @constructor
         */
        this.MyEdge = function (i, j, x1, x2) {
            this.i = i;
            this.j = j;
            this.x1 = x1;
            this.x2 = x2;
        };
        var U = [], B = [];
    //U skeleton nodes :vector3 array ,
    // B skeleton edges :edge2 array  .:
    //contraction: an object of class Contraction
    // connection; an object of class ConnectivitySurgery
    // regfinement

};

THREE.SkeletonGeometry.prototype = {
    constructor: THREE.SkeletonGeometry,
    /**
     * find all the Face3s adjacent to i, and return a array storing those Face3s of all the other vertexes
     * @param i
     * @param size
     * @param F
     * @returns {Array}
     */
    findAllFace3: function (i, F) {
        "use strict";
        var faces = [],
            count,
            count2 = 0;
        for (count = 0; count < F.length; count += 1) {
            if (F[count].a === i || F[count].b === i || F[count].c === i) {
                faces[count2] = F[count];
                count2 += 1;
            }
        }
        return faces;
    },
    /**
     * findAllEdges: find all the edges related to i and return an array of MyEdge
     * @param i
     * @param faces
     * @returns {Array}
     */
    findAllEdges: function (i, faces) {
            "use strict";
            var count, count2 = 0,
                edges = [],
                flags = [],
                x, y,
                undiscovered = "UNDISCOVERED", discovered = "DISCOVERED",
                isContained = function (i, face) {
                    switch(i) {
                        case face.a: return 0;
                        case face.b: return 1;
                        case face.c: return 2;
                        default : return null;
                    }
                },//return null if  face doesn't contain i, or return the 0, 1 or 2 (a,b,c)
                findTheLastVertex = function (i ,j, x) {
                    var positionI = null,
                        positionJ = null;
                    for (count = 0; count < faces.length; count += 1) {
                            positionI = isContained(i, faces[count]);
                            positionJ = isContained(j, faces[count]);
                            if( (positionI !== null) && (positionJ !== null)) {
                                flags[count] = discovered;
                                switch(positionI + positionJ) {
                                    case 1: if (faces[count].c !== x) {
                                        return faces[count].c;}
                                        break;
                                    case 2: if (faces[count].b !== x) {
                                        return faces[count].b;}
                                        break;
                                    case 3: if (faces[count].a !== x) {
                                        return faces[count].a;}
                                        break;
                                    default: return null;
                                }
                            }
                    }
                };//return the last vertexIndex.  return null if error
            for (count = 0; count < faces.length; count += 1) {
                flags[count] = undiscovered;
            }
            for (count = 0; count < faces.length; count += 1) {
                if(flags[count] === undiscovered) {
                    switch (i) {
                        case faces[count].a:
                            x = faces[count].b;
                            y = faces[count].c;
                            break;
                        case faces[count].b:
                            x = faces[count].a;
                            y = faces[count].c;
                            break;
                        case faces[count].c:
                            x = faces[count].a;
                            y = faces[count].b;
                            break;
                    } //x and y are the other 2 coplanar vertexes
                    flags[count] = discovered;
                    edges[count2] = new this.MyEdge(i, x, y, findTheLastVertex(i, x, y));
                    count2 += 1;
                    edges[count2] = new this.MyEdge(i, y, x, findTheLastVertex(i, y, x));
                    count2 += 1;
                }
            }
            return edges;
        },
    generateSkeleton: function () {
        var contraction = new THREE.SkeletonGeometry.Contraction(this.geo.vertices, this.geo.faces);
        contraction.contract();
    },
    /**
     * return true if A equals B approximately
     * @param A
     * @param B
     * @returns {boolean}
     */
    compareMatrices: function (A, B) {
        var x = A.length,
            y = A[0].length,
            count, count2;
        if (A.length !== B.length || A[0].length !== B[0].length) {
            return false;
        }
        for (count = 0; count < x; count += 1) {
            for (count2 = 0; count2 < y; count2 += 1) {
                if (Math.round(A[count][count2]) !== Math.round(B[count][count2])) {
                    return false;
                }
            }
        }
        return true;

    }
};

THREE.SkeletonGeometry.LaplaceMatrix = function(V, F) {//  LaplaceMatrix can be  used to compute the Laplace matrix L according to current vertices vector V
    this.V = V.concat();
    this.F = F.concat();
    this.Vsize = V.length;
    this.Fsize = F.length;
    this.LMatrix = [];

};

THREE.SkeletonGeometry.LaplaceMatrix.prototype = {
    constructor: THREE.SkeletonGeometry.LaplaceMatrix,
    /**
     * set 0 for each member of the matrix first
     */
    initLMatrix: function () {

        "use strict";
        var count, count2;
        for (count = 0; count < this.Vsize; count += 1) {
            this.LMatrix[count] = [];
            for (count2 = 0; count2 < this.Vsize; count2 += 1) {
                this.LMatrix[count][count2] = 0;
            }
        }
        return this.LMatrix;

    },
    /**
     * compute Wij1
     * @param myEdge
     * @returns {number}
     */
    getWij1: function (myEdge) {
        var ij, ix1, ix2, jx1, jx2, cosX1, cosX2, cotX1, cotX2;
        ij = new THREE.Line3(this.V[myEdge.i],this.V[myEdge.j]).distance();
        ix1 = new THREE.Line3(this.V[myEdge.i],this.V[myEdge.x1]).distance();
        ix2 = new THREE.Line3(this.V[myEdge.i],this.V[myEdge.x2]).distance();
        jx1 = new THREE.Line3(this.V[myEdge.j],this.V[myEdge.x1]).distance();
        jx2 = new THREE.Line3(this.V[myEdge.j],this.V[myEdge.x2]).distance();
        cosX1 = (ix1 * ix1 + jx1 * jx1 - ij * ij) / (2 * ix1 * jx1);
        cosX2 = (ix2 * ix2 + jx2 * jx2 - ij * ij) / (2 * ix2 * jx2);
        cotX1 = cosX1 / Math.sqrt(1 - cosX1 * cosX1);
        cotX2 = cosX2 / Math.sqrt(1 - cosX2 * cosX2);
        return cotX1 + cotX2;
    },

    /**
     * set all the wij2 for this.LMatrix
     */
    setWij2: function (){
        var i, j, wij2;
        for (i = 0; i < this.Vsize; i += 1) {
            wij2 = 0;
            for (j = 0; j < this.Vsize; j += 1) {
                if (i !== j) {
                    wij2 -= this.LMatrix[i][j];
                }
            }
            this.LMatrix[i][i] = wij2;
        }
    },

    /**
     * generate Laplace matrix according to current V
     */
    generateL: function () {
        var count, count2,
            edges,faces,
            wij1,
            vertex1,
            vertex2,
            SG;
        SG = new THREE.SkeletonGeometry();
        this.initLMatrix();
        for (count = 0; count < this.Vsize; count += 1) {
            faces = SG.findAllFace3(count, this.F);
            edges = SG.findAllEdges(count, faces);
            vertex1 = count;
            for (count2 = 0; count2 < edges.length; count2 += 1) {
                vertex2 = edges[count2].j;
                wij1 = this.getWij1(edges[count2]);
                this.LMatrix[vertex1][vertex2] = wij1;
                this.LMatrix[vertex2][vertex1] = wij1;
            }
        }
        this.setWij2();
        return this.LMatrix;
    },

    setVF: function (V, F) {
        this.V = V;
        this.F = F;
    }
};

THREE.SkeletonGeometry.Contraction = function(V, F) {// Contraction can be used to  contract the original geometry
    this.V = V.concat();
    this.F = F.concat();
    this.WH = [];
    this.WH0 = null;
    this.AO = [];
    this.A = [];
    this.SL = 2.0;
    this.L = null;
    this.WL = null;
};

THREE.SkeletonGeometry.Contraction.prototype = {
    construcor: THREE.SkeletonGeometry.Contraction,

    contract: function () {
        var counter = 0;
        this.L = new THREE.SkeletonGeometry.LaplaceMatrix(this.V, this.F);
        this.initWeight();
        while(this.computeS() > 1 && counter !== 1) {
            this.L.setVF(this.V, this.F);
            this.L.generateL();
            this.setNextV();
            this.updateWeight();
            counter += 1;
        }
        console.log("counter = ", counter);
    },

    /**
     * return the total superficial area
     * @param V
     * @param F
     * @returns {number}
     */
    computeS: function () {
        var S = 0,
            count;
        for (count = 0;count < this.F.length; count += 1) {
            S += this.computeSingleArea(this.F[count]);
        }
        return S;
    },
    computeSingleArea: function (face) {
        var tri = new THREE.Triangle(this.V[face.a], this.V[face.b], this.V[face.c]);
        return tri.area();

    },
    getOneRingArea: function (i) {
        var SG = new THREE.SkeletonGeometry(),
            faces,
            OFA = 0,
            count;
        faces = SG.findAllFace3(i, this.F);
        for (count = 0; count < faces.length; count += 1) {
            OFA += this.computeSingleArea(faces[count]);
        }
        return OFA;

    },
    initWeight: function () {
        var count;
        this.WH0 = 1;
        for (count = 0; count < this.V.length; count += 1) { //initiate WH[]
            this.WH[count] = this.WH0;
        }
        for (count = 0; count < this.V.length; count += 1) { //initiate A0[]
            this.AO[count] = this.getOneRingArea(count);
        }
        this.A = this.AO.concat();
        this.WL = 0.001 * Math.sqrt(this.computeS() / this.F.length); // initiate WL
    },
    updateWeight: function () {
        var count;
        this.WL *= this.SL; //update WL
        for (count = 0; count < this.V.length; count += 1) { //update A[]
            this.A[count] = this.getOneRingArea(count);
        }

        for (count = 0; count < this.V.length; count += 1){ //update WH[]
            this.WH[count] = this.WH0 * Math.sqrt(this.AO[count] / this.A[count]);
        }

//        console.info("after updateWeight()\n", this);
    },
    setNextV: function () {
        var A, B, A1, B1, transA1, A2, B2, inA2, count, VV;
        A = this.setMatrixA();
        B = this.setMatrixB();
        A1 = new Matrix();
        B1 = new Matrix();
        A1.setElements(A);
        B1.setElements(B);
        transA1 = A1.transpose();
//        console.log(A1, B1);
        A2 = transA1.multiply(A1);
        B2 = transA1.multiply(B1);
//        console.log("A2 = ", A2);


        inA2 = new Matrix.create(numeric.inv(A2.elements));
        console.log("inA2 = ", inA2);
        VV = inA2.multiply(B2).elements;
//        A2.print();
//        inA2.print();
//        B2.print();
        console.log("VV = ", VV);
        if(inA2 === null) {
            console.alert("Matrix A2 doesn't have an inverse matrix!");
        }
        for (count = 0; count < this.V.length; count += 1) {
            this.V[count] = new THREE.Vector3(VV[count][0], VV[count][1], VV[count][2]);
        }
        return VV;
    },

    setMatrixA: function () {
        var A = [this.V.length + 1],
            count, count2;
        for (count = 0; count < this.V.length; count += 1) {
            A[count] = [this.V.length];
            for (count2 = 0; count2 < this.V.length; count2 += 1) {
                A[count][count2] = this.L.LMatrix[count][count2] * this.WL;
            }
        }
        A[this.V.length] = [this.V.length];
        for (count = 0; count < this.V.length; count += 1 ) {
            A[this.V.length][count] = this.WH[count];
        }
        return A;
    },
    setMatrixB: function () {
        var B = [this.V.length + 1],
            count, count2,
            X, Y, Z, V;
        for (count = 0; count < this.V.length; count += 1) {
            B[count] = [3];
            for (count2 = 0; count2 < 3; count2 += 1) {
                B[count][count2] = 0;
            }
        }
        X = new Matrix();
        X.setElements(this.WH);
        V = [this.V.length];
        for (count = 0; count < this.V.length; count += 1) {
            V[count] = [3];
            V[count][0] = this.V[count].x;
            V[count][1] = this.V[count].y;
            V[count][2] = this.V[count].z;
        }
        Y = new Matrix();
        Y.setElements(V);
        X = X.transpose();
        Z = X.multiply(Y);
        B[this.V.length] = [];
        for (count2 = 0; count2 < 3; count2 += 1) {
            B[this.V.length][count2] = Z.elements[0][count2];
        }
        return B;
    }
};