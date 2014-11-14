/**
 * Created by mimicry on 2014/11/10.
 */

module("SkeletonGeometry");

//test("initLMatrix()" , function () {
//    var cube = new THREE.BoxGeometry(1,1,1),
//        laplace = new THREE.SkeletonGeometry.LaplaceMatrix(cube.vertices, cube.faces),
//        SG =  new THREE.SkeletonGeometry(cube);
//        console.info(cube);
////    console.info(laplace);
//    var emptyL = laplace.initLMatrix();
//
//    var emptyMatrix = [],
//        i, j;
//
//    for (i = 0; i < 8; i += 1) {
//        emptyMatrix[i] = [];
//        for (j = 0; j < 8; j += 1) {
//            emptyMatrix[i][j] = 0;
//        }
//    }
////    console.info(emptyMatrix);
//    ok(SG.compareMatrices(emptyL, emptyMatrix), "Passed!");
//});
//test("findAllFace3()", function () {
//    var cube = new THREE.BoxGeometry(1,1,1),
//        LM = new THREE.SkeletonGeometry.LaplaceMatrix(cube.vertices, cube.faces),
//        SG = new THREE.SkeletonGeometry(cube),
//        correctFaces = [LM.F[0], LM.F[5], LM.F[8], LM.F[9]],
//        faces,
//        isEqual = true,
//        count;
//    faces = SG.findAllFace3(0, LM.F);
////    console.info(correctFaces);
////    console.info(faces);
//    for(count = 0; count < faces.length; count += 1) {
//        if (faces[count] !== correctFaces[count]) {
//            isEqual = false;
//            break;
//        }
//    }
//    ok(isEqual, "Passed!");
//});
//test("findAllEdges()", function () {
//    var cube = new THREE.BoxGeometry(1,1,1),
//        LM = new THREE.SkeletonGeometry.LaplaceMatrix(cube.vertices, cube.faces),
//        SG = new THREE.SkeletonGeometry(cube),
//        faces, edges,
//        isEqual = true,
//        count;
//    faces = SG.findAllFace3(0, LM.F);
//    console.info(faces);
//    edges = SG.findAllEdges(0, faces);
//    console.info(edges);
//    ok(isEqual, "Passed!");
//});
//test("generateL()", function () {
//    var cube, LM, SG, correctMatrix, L;
//    cube = new THREE.BoxGeometry(1,1,1);
//    LM = new THREE.SkeletonGeometry.LaplaceMatrix(cube.vertices, cube.faces);
//    SG = new THREE.SkeletonGeometry(cube);
//    correctMatrix = [[-6, 2, 2, 0, 0, 2, 0, 0],
//                    [2, -6, 0, 2, 2, 0, 0, 0],
//                    [2, 0, -6, 2, 0, 0, 0, 2],
//                    [0, 2, 2, -6, 0, 0, 2, 0],
//                    [0, 2, 0, 0, -6, 2, 2, 0],
//                    [2, 0, 0, 0, 2, -6, 0, 2],
//                    [0, 0, 0, 2, 2, 0, -6, 2],
//                    [0, 0, 2, 0, 0, 2, 2, -6]];
//    L = LM.generateL();
//    console.info(L);
//    ok(SG.compareMatrices(correctMatrix, L), "Passed!");
//});

test("initWeight()", function () {
    var cube, C;
    cube = new THREE.BoxGeometry(1,1,1);
    C = new THREE.SkeletonGeometry.Contraction(cube.vertices, cube.faces);
    C.initWeight();
    ok(true , "Passed!");
});

test("computeS()", function () {
    var cube, C;
    cube = new THREE.BoxGeometry(1,1,1);
    C = new THREE.SkeletonGeometry.Contraction(cube.vertices, cube.faces);
    ok(C.computeS() === 6 , "Passed!");
});

test("setNextV()", function () {
    var cube, C, L, count, correct = [],
        count2, SG;
    cube = new THREE.BoxGeometry(1,1,1);
    SG = new THREE.SkeletonGeometry(cube);
    for (count = 0; count < 8; count += 1) {
        cube.vertices[count].z += 1;
        correct[count] = [];
        correct[count][0] = 0;
        correct[count][1] = 0;
        correct[count][2] = 1;
    }
    C = new THREE.SkeletonGeometry.Contraction(cube.vertices, cube.faces);
    L = new THREE.SkeletonGeometry.LaplaceMatrix(cube.vertices, cube.faces);
    C.L = L;
    C.L.generateL();
    C.initWeight();

    ok(SG.compareMatrices(C.setNextV(), correct), "Passed!");
});

test("updateWeight()", function () {
    var cube, C, L, count, correct, count2;
    cube = new THREE.BoxGeometry(2, 2, 2);
    C = new THREE.SkeletonGeometry.Contraction(cube.vertices, cube.faces);
    L = new THREE.SkeletonGeometry.LaplaceMatrix(cube.vertices, cube.faces);
    C.L = L;
    C.L.generateL();
    C.initWeight();
    C.setNextV();
    C.updateWeight();
    console.info(C);
    ok(true , "Passed!");
});

test("contract()", function () {
    var cube, C, L, count, correct = [], count2;
    cube = new THREE.BoxGeometry(2, 2, 2);
    for (count = 0; count < 8; count += 1) {
        cube.vertices[count].z += 1;
        correct[count] = [];
        correct[count][0] = 0;
        correct[count][1] = 0;
        correct[count][2] = 1;
    }
    C = new THREE.SkeletonGeometry.Contraction(cube.vertices, cube.faces);
    C.contract();
    console.info(C);
    ok(true , "Passed!");
});