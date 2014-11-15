/**
 * Created by mimicry on 2014/11/8.
 */

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
    this.WH = new Matrix();
    this.WH0 = new Matrix();
    this.AO = [];
    this.A = [];
    this.SL = 2.0;
    this.L = null;
    this.WL = new Matrix();
};

THREE.SkeletonGeometry.Contraction.prototype = {
    constructor: THREE.SkeletonGeometry.Contraction,

    contract: function () {
        var counter = 0;
        this.L = new THREE.SkeletonGeometry.LaplaceMatrix(this.V, this.F);
        this.initWeight();
        while(this.computeS() > 100) {
            this.L.setVF(this.V, this.F);
            this.L.generateL();
//            console.info("Matrix L:\n");
//            Matrix.create(this.L.LMatrix).print();
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
        var count, array = [], wl;


        for (count = 0; count < this.V.length; count += 1) { //initiate WH0[]
            array[count] = 1.0;
        }
        this.WH0 = Matrix.Diagonal(array);
        this.WH = this.WH0.dup();

        for (count = 0; count < this.V.length; count += 1) { //initiate A0[]
            this.AO[count] = this.getOneRingArea(count);
        }
        this.A = this.AO.concat();

        wl = 0.001 * Math.sqrt(this.computeS() / this.F.length); // initiate WL
        array = [];
        for (count = 0; count < this.V.length; count += 1) { //initiate WH0[]
            array[count] = wl;
        }
        this.WL= Matrix.Diagonal(array);
    },
    updateWeight: function () {
        var count, array = [];
        this.WL = this.WL.multiply(this.SL); //update WL
        for (count = 0; count < this.V.length; count += 1) { //update A[]
            this.A[count] = this.getOneRingArea(count);
        }

        for (count = 0; count < this.V.length; count += 1){ //update WH[]
            array[count] = this.WH0.e(1, 1) * Math.sqrt(this.AO[count] / this.A[count]);
        }
        this.WH = Matrix.Diagonal(array);

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
//        A1.print();
//        B1.print();
        transA1 = A1.transpose();
//        console.log(A1, B1);
        A2 = transA1.multiply(A1);
        B2 = transA1.multiply(B1);
//        console.log("A2 = ", A2);


        inA2 = new Matrix.create(numeric.inv(A2.elements));
//        console.log("inA2 = ", inA2);
        VV = inA2.multiply(B2).elements;
//        A2.print();
//        inA2.print();
//        B2.print();
//        console.log("VV = ", VV);
        if(inA2 === null) {
            console.alert("Matrix A2 doesn't have an inverse matrix!");
        }
        for (count = 0; count < this.V.length; count += 1) {
            this.V[count] = new THREE.Vector3(VV[count][0], VV[count][1], VV[count][2]);
        }
        return VV;
    },

    setMatrixA: function () {
        var A = [this.V.length * 2],
            count, count2, wll, l;
        l = new Matrix.create(this.L.LMatrix);

        wll = this.WL.multiply(l);
        for (count = 0; count < this.V.length; count += 1) {
            A[count] = [this.V.length];
            for (count2 = 0; count2 < this.V.length; count2 += 1) {
//                console.info("e( ", count , count2, ") = ", wll.e(count + 1, count2 + 1));
                A[count][count2] = wll.e(count + 1, count2 + 1);
            }
        }

        for (count = this.V.length; count < 2 * this.V.length; count += 1 ) {
            A[count] = [this.V.length];
            for (count2 = 0; count2 < this.V.length; count2 += 1) {
                A[count][count2] = this.WH.e(count - this.V.length + 1, count2 + 1);
            }
        }
        return A;
    },
    setMatrixB: function () {
        var B = [this.V.length * 2],
            count, count2,
            X, Y, Z, V;
        for (count = 0; count < this.V.length; count += 1) {
            B[count] = [3];
            for (count2 = 0; count2 < 3; count2 += 1) {
                B[count][count2] = 0;
            }
        }
        X = this.WH.dup();
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

        for (count = this.V.length; count < this.V.length * 2; count += 1) {
            B[count] = [3];
            for (count2 = 0; count2 < 3; count2 += 1) {
                B[count][count2] = Z.e(count - this.V.length + 1, count2 + 1);
            }
        }
        return B;
    }
};