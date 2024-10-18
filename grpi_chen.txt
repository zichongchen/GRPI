var polys = 
    /* color: #d63000 */
    /* displayProperties: [
      {
        "type": "rectangle"
      },
      {
        "type": "rectangle"
      }
    ] */
    ee.Geometry.MultiPolygon(
      [[[[-122, 39],
         [-122, 37],
         [-120, 37],
         [-120,39]]],
         [[[-122, 37],
         [-122, 35],
         [-120, 35],
         [-120,37]]],
         [[[-120, 39],
         [-120, 37],
         [-118, 37],
         [-118,39]]],
         [[[-120, 37],
         [-120, 35],
         [-118, 35],
         [-118, 37]]]],null,false);



polys.coordinates().size().evaluate(function (count) {
  for (var i= 0; i < count; i++) { // i is index of the polygon
    // Extract the coordinates for the polygon
    var polygon = ee.Geometry.Polygon(
              ee.List(polys.coordinates().get(i)) 
            )
            // Do work with the simple polygon
            doTropwet(polygon, i)
          }
        })


function doTropwet(geometry, poly_ind) {
  // Most of your code goes here
  Map.centerObject(geometry) 
  
  // set year and month range
  var first_year = 2022;
  var last_year = 2022; 
  var first_month =9;
  var last_month = 9;

  // scale factors
  var parScale = 16
  var tileScale = 16
  var resolution = 120
  
  // switch in case specific missions want to be used
  var expressionCK9 = 1
  var expressionCK8 = 1
  var expressionCK7 = 1
  var expressionCK5 = 1

  // region geometry used interchangeably
  var region = geometry 
  var geometry = geometry
  
  // add geometry to map
  Map.addLayer(geometry,{},"Total Processing Area Geometry")
  
  //// calculate Slope and fetch HAND \\\\
  
  // var hand = ee.Image('users/gena/GlobalHAND/90m-global/hand-1000')
  var hand = ee.Image('users/gena/GlobalHAND/30m/hand-1000')

  var dataset = ee.Image('CGIAR/SRTM90_V4');
  var elevation = dataset.select('elevation');
  var elevation = elevation.clip(geometry)
  var elevation = elevation.unmask(0)
  var slope = ee.Terrain.slope(elevation);
  var slope = slope.rename(['slope'])  
          
  var slopeMask = slope.lte(2.5)

  var handMask = hand.lt(15)  

  var handSlopeMask = ee.Image(slopeMask.add(handMask)).eq(2)

////   
//// functions for implementing cloud/QA masks
//// 

// Landsat 8
  
    var cloudMaskL8 = function(image) {
    var cloudShadowBitMask = (1 << 3);
    var cloudsBitMask = (1 << 5);
    var qa = image.select('QA_PIXEL').clip(geometry);
    var ra = image.select('QA_RADSAT')
    /// Check that the cloud bit is off.
    // See https://landsat.usgs.gov/collectionqualityband
    var mask = qa.bitwiseAnd(cloudShadowBitMask).eq(0)
        .and(qa.bitwiseAnd(cloudsBitMask).eq(0))
        .or(ra.bitwiseAnd(1<<1))
    
    var maskedImage = image.updateMask(mask);
    
    var bandNames = ['SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B6', 'SR_B7'];
    var maskedImageSelectedBands = maskedImage.select(bandNames);
    maskedImageSelectedBands = maskedImageSelectedBands.rename('b','g','r','nir','swir1','swir2')
    maskedImageSelectedBands = maskedImageSelectedBands.short()
    return maskedImageSelectedBands.clip(geometry);
  };

// Landsat 9
  var cloudMaskL9 = function(image) {
    var cloudShadowBitMask = (1 << 3);
    var cloudsBitMask = (1 << 5);
    var qa = image.select('QA_PIXEL').clip(geometry);
    var ra = image.select('QA_RADSAT')
    /// Check that the cloud bit is off.
    // See https://landsat.usgs.gov/collectionqualityband
    var mask = qa.bitwiseAnd(cloudShadowBitMask).eq(0)
        .and(qa.bitwiseAnd(cloudsBitMask).eq(0))
        .or(ra.bitwiseAnd(1<<1))
    
    var maskedImage = image.updateMask(mask);
    
    var bandNames = ['SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B6', 'SR_B7'];
    var maskedImageSelectedBands = maskedImage.select(bandNames);
    maskedImageSelectedBands = maskedImageSelectedBands.rename('b','g','r','nir','swir1','swir2')
    maskedImageSelectedBands = maskedImageSelectedBands.short()
    return maskedImageSelectedBands.clip(geometry);
  };

// Landsat 4,5,7  
  var cloudMaskL457 = function(image) {
    var qa = image.select('QA_PIXEL').clip(geometry);
    var ra = image.select('QA_RADSAT')
    /// Check that the cloud bit is off.
    // See https://landsat.usgs.gov/collectionqualityband
    var qaMask = qa.bitwiseAnd(parseInt('11111',2)).eq(0);
    var satMask = ra.eq(0);
    
    // var maskedImage = image.updateMask(mask);
    var maskedImage = image.updateMask(qaMask).updateMask(satMask);
    
    var bandNames = ['SR_B1', 'SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B7'];
    var maskedImageSelectedBands = maskedImage.select(bandNames);
    maskedImageSelectedBands = maskedImageSelectedBands.rename('b','g','r','nir','swir1','swir2')

    
    
    maskedImageSelectedBands = maskedImageSelectedBands.short()
    return maskedImageSelectedBands.clip(geometry);
  };
  
  // renaming bands
  
  var genFalseND = function(){
  
    var b = ee.Image(ee.Number(0)).clip(geometry).cast({"constant":"short"}).rename("b")
    var b = b.addBands(ee.Image(ee.Number(0)).clip(geometry).cast({"constant":"short"}).rename("g"))
    var b = b.addBands(ee.Image(ee.Number(0)).clip(geometry).cast({"constant":"short"}).rename("r"))
    var b = b.addBands(ee.Image(ee.Number(0)).clip(geometry).cast({"constant":"short"}).rename("nir"))
    var b = b.addBands(ee.Image(ee.Number(0)).clip(geometry).cast({"constant":"short"}).rename("swir1"))
    var b = b.addBands(ee.Image(ee.Number(0)).clip(geometry).cast({"constant":"short"}).rename("swir2"))
    var b = b.updateMask(b)
    return ee.Image(b)
  
}

  print("Get Landsat Data..."+poly_ind)

  //// Functions to Fetch LS8 and LS 9 imagery and implement cloud mask \\\\

   var getLS9 = function(){

    var expression = ee.Number(first_month).gt(ee.Number(last_month))
    var firstHalfLS9 = ee.Algorithms.If(expression, ee.ImageCollection("LANDSAT/LC09/C02/T1_L2").filterBounds(geometry).filter(ee.Filter.calendarRange(first_year,first_year,'year')).filter(ee.Filter.calendarRange(first_month,12,'month')),ee.ImageCollection('LANDSAT/LC09/C02/T1_L2').filterDate('1910-01-01', '1910-12-31'))
    var secondHalfLS9 = ee.Algorithms.If(expression, ee.ImageCollection("LANDSAT/LC09/C02/T1_L2").filterBounds(geometry).filter(ee.Filter.calendarRange(last_year,last_year,'year')).filter(ee.Filter.calendarRange(1,last_month,'month')), ee.ImageCollection("LANDSAT/LC09/C02/T1_L2").filterBounds(geometry).filter(ee.Filter.calendarRange(first_year,last_year,'year')).filter(ee.Filter.calendarRange(first_month,last_month,'month')))
    var dateSelectLS9 = ee.ImageCollection(ee.ImageCollection(firstHalfLS9).merge(ee.ImageCollection(secondHalfLS9)))
    // print(dateSelectLS8.size())
    var conditional = dateSelectLS9.size().gt(0)
    var ls9CloudMasked = ee.Algorithms.If(conditional,dateSelectLS9.map(cloudMaskL9),ee.ImageCollection.fromImages([genFalseND()]))

    // var ls8CloudMasked = ee.ImageCollection(dateSelectLS8.map(cloudMaskL8))
    return ee.ImageCollection(ls9CloudMasked) 
  }
  var ls9CloudMasked = ee.ImageCollection(ee.Algorithms.If(expressionCK9,getLS9(),ee.ImageCollection('LANDSAT/LC09/C02/T1_L2').filterDate('1910-01-01', '1910-12-31')))

  
  var getLS8 = function(){

    var expression = ee.Number(first_month).gt(ee.Number(last_month))
    var firstHalfLS8 = ee.Algorithms.If(expression, ee.ImageCollection("LANDSAT/LC08/C02/T1_L2").filterBounds(geometry).filter(ee.Filter.calendarRange(first_year,first_year,'year')).filter(ee.Filter.calendarRange(first_month,12,'month')),ee.ImageCollection('LANDSAT/LC08/C02/T1_L2').filterDate('1910-01-01', '1910-12-31'))
    var secondHalfLS8 = ee.Algorithms.If(expression, ee.ImageCollection("LANDSAT/LC08/C02/T1_L2").filterBounds(geometry).filter(ee.Filter.calendarRange(last_year,last_year,'year')).filter(ee.Filter.calendarRange(1,last_month,'month')), ee.ImageCollection("LANDSAT/LC08/C02/T1_L2").filterBounds(geometry).filter(ee.Filter.calendarRange(first_year,last_year,'year')).filter(ee.Filter.calendarRange(first_month,last_month,'month')))
    var dateSelectLS8 = ee.ImageCollection(ee.ImageCollection(firstHalfLS8).merge(ee.ImageCollection(secondHalfLS8)))
    // print(dateSelectLS8.size())
    var conditional = dateSelectLS8.size().gt(0)
    var ls8CloudMasked = ee.Algorithms.If(conditional,dateSelectLS8.map(cloudMaskL8),ee.ImageCollection.fromImages([genFalseND()]))

    // var ls8CloudMasked = ee.ImageCollection(dateSelectLS8.map(cloudMaskL8))
    return ee.ImageCollection(ls8CloudMasked) 
  }
  var ls8CloudMasked = ee.ImageCollection(ee.Algorithms.If(expressionCK8,getLS8(),ee.ImageCollection('LANDSAT/LC08/C02/T1_L2').filterDate('1910-01-01', '1910-12-31')))

  
 //// Functions to Fetch LS7 imagery and implement cloud mask \\\\

  var getLS7 = function(){

    var expression = ee.Number(first_month).gt(ee.Number(last_month))
    var firstHalfLS7 = ee.Algorithms.If(expression, ee.ImageCollection("LANDSAT/LE07/C02/T1_L2").filterBounds(geometry).filter(ee.Filter.calendarRange(first_year,first_year,'year')).filter(ee.Filter.calendarRange(first_month,12,'month')),ee.ImageCollection('LANDSAT/LE07/C02/T1_L2').filterDate('1910-01-01', '1910-12-31'))
    var secondHalfLS7 = ee.Algorithms.If(expression, ee.ImageCollection("LANDSAT/LE07/C02/T1_L2").filterBounds(geometry).filter(ee.Filter.calendarRange(last_year,last_year,'year')).filter(ee.Filter.calendarRange(1,last_month,'month')), ee.ImageCollection("LANDSAT/LE07/C02/T1_L2").filterBounds(geometry).filter(ee.Filter.calendarRange(first_year,last_year,'year')).filter(ee.Filter.calendarRange(first_month,last_month,'month')))
    var dateSelectLS7 = ee.ImageCollection(firstHalfLS7).merge(secondHalfLS7)
    // print(dateSelectLS7.size())
    var conditional = dateSelectLS7.size().gt(0)
    var ls7CloudMasked = ee.Algorithms.If(conditional,dateSelectLS7.map(cloudMaskL457),ee.ImageCollection.fromImages([genFalseND()]))
    
    // var ls7CloudMasked = dateSelectLS7.map(cloudMaskL457)
    return ee.ImageCollection(ls7CloudMasked)
  }
  var ls7CloudMasked = ee.ImageCollection(ee.Algorithms.If(expressionCK7, getLS7(),ee.ImageCollection('LANDSAT/LE07/C02/T1_L2').filterDate('1910-01-01', '1910-12-31')))

//// 
// Function to fill gaps due to SLC error using focal mean 
// NOTE: this will also fill gaps in LS7 composite due to cloud cover
  function fillGap(image){
  return image.focal_mean(1, 'square', 'pixels', 8).blend(image);
  }

  
  //// Function to Fetch and Correct LS5
  var getLS5 = function(){

    var expression = ee.Number(first_month).gt(ee.Number(last_month))
    var firstHalfLS5 = ee.Algorithms.If(expression, ee.ImageCollection("LANDSAT/LT05/C02/T1_L2").filterBounds(geometry).filter(ee.Filter.calendarRange(first_year,first_year,'year')).filter(ee.Filter.calendarRange(first_month,12,'month')),ee.ImageCollection('LANDSAT/LT05/C02/T1_L2').filterDate('1910-01-01', '1910-12-31'))
    var secondHalfLS5 = ee.Algorithms.If(expression, ee.ImageCollection("LANDSAT/LT05/C02/T1_L2").filterBounds(geometry).filter(ee.Filter.calendarRange(last_year,last_year,'year')).filter(ee.Filter.calendarRange(1,last_month,'month')), ee.ImageCollection("LANDSAT/LT05/C02/T1_L2").filterBounds(geometry).filter(ee.Filter.calendarRange(first_year,last_year,'year')).filter(ee.Filter.calendarRange(first_month,last_month,'month')))
    
    // print(firstHalfLS5)
    
    var dateSelectLS5 = ee.ImageCollection(firstHalfLS5).merge(secondHalfLS5)
    // print(dateSelectLS5.size())
    var conditional = dateSelectLS5.size().gt(0)
    var ls5CloudMasked = ee.Algorithms.If(conditional,dateSelectLS5.map(cloudMaskL457),ee.ImageCollection.fromImages([genFalseND()]))
    
    var ls5CloudMasked = dateSelectLS5.map(cloudMaskL457)
    return ee.ImageCollection(ls5CloudMasked) 
  }
  var ls5CloudMasked = ee.ImageCollection(ee.Algorithms.If(expressionCK5,getLS5(),ee.ImageCollection('LANDSAT/LT05/C02/T1_L2').filterDate('1910-01-01', '1910-12-31')))

  
  //// Generate Composites \\\\
  
  // .merge(ee.ImageCollection(ls7CloudMasked))

  print("Gen Comp")
  
  print('LS5: ',ls5CloudMasked)
  print('LS7: ',ls7CloudMasked)
  print('LS8: ',ls8CloudMasked)
  print('LS9: ',ls9CloudMasked)
  // print('LS9: ',ls8CloudMasked)
  
  var median_compositeCol = ls5CloudMasked.merge(ls8CloudMasked).merge(ls7CloudMasked).merge(ls9CloudMasked)
  
// 
// composite original, non-gap filled LS scenes for comparison
// 
 var median_composite_orig = median_compositeCol.reduce(ee.Reducer.median(),parScale) 
 //Map.addLayer(median_composite_orig, {bands: 'nir_median,swir1_median,b_median', min:5000.0,  max: 24000, gamma: 1.4}, 'Landsat composite original');
// 
// 
//
  
  
  
  var median_compositeCol = median_compositeCol.map(fillGap)
  
  // print(median_compositeCol)
  
  var median_composite = median_compositeCol.reduce(ee.Reducer.median(),parScale)
  
  var median_composite = median_composite.select(
    ['b_median', 'g_median', 'r_median','nir_median','swir1_median','swir2_median'], // old names
    ['b', 'g', 'r','nir','swir1','swir2'])
    
  var median_composite_scaled = median_composite.multiply(0.0000275).add(-0.2);
  var median_composite_scaled = median_composite_scaled.multiply(10000);

  // view landsat composite
  //Map.addLayer(median_composite_scaled, {bands: 'nir,swir1,b', min:0,  max: 3000, gamma: 1.4}, 'Landsat composite scaled');
  
  
  ////
  //// automatically generate endmemebers
  ////
  var lsuImage = median_composite_scaled;
  
  // calc NDWI
  var ndiiComp = lsuImage.normalizedDifference(['nir','swir1'])
  
  // function for running ruleset to extract endmembers
  // for water, veg and bare soil
  var genEndmembers = function(lsuImage){
    // NDVI
    var ndviComp = lsuImage.normalizedDifference(['nir', 'r']);
    //Map.addLayer(ndviComp, {}, 'NDVI');
    // mNDWI
    var ndwiComp = lsuImage.normalizedDifference(['g', 'swir1']);
    //Map.addLayer(ndwiComp, {}, 'NDWI');
    
    // Run conditional statements on NDVI to train pure pixel collection.
    // Create binary layers using logical operations:
    var waterThresh = ndwiComp.gt(0.5);
    var gvThresh = ndviComp.gt(0.7);
    //var bsThresh = ndwiComp.gt(-0.285).and(ndwiComp.lt(-0.26));
    var bsThresh2 = ndwiComp.gt(-0.285).and(ndwiComp.lt(-0.26).and(ndviComp.gt(0.16).and(ndviComp.lt(0.17)))); 
    
    // Apply masks to composite image to create three individual images:
    var lsuImageWater = lsuImage.updateMask(waterThresh);
    //Map.addLayer(lsuImageWater, {}, 'pure_water',false);
    var lsuImageGV = lsuImage.updateMask(gvThresh);
    //Map.addLayer(lsuImageGV, {}, 'pure_veg',false);
    var lsuImageBS = lsuImage.updateMask(bsThresh2);
    //Map.addLayer(lsuImageBS, {}, 'pure_sand',false);
    
    // Find mean reflectance across the 6 bands for each masked image:
    
    var waterMean = lsuImageWater.reduceRegion({
      reducer: ee.Reducer.mean(),
      geometry: region,
      scale: resolution,
      tileScale: tileScale,
      maxPixels: 3e10
    });
    
    var waterMean = ee.List([waterMean.get('b'),waterMean.get('g'),waterMean.get('r'),waterMean.get('nir'),waterMean.get('swir1'),waterMean.get('swir2')]);
    
    var vegMean = lsuImageGV.reduceRegion({
      reducer: ee.Reducer.mean(),
      geometry: region,
      scale: resolution,
      tileScale: tileScale,
      maxPixels: 3e10
    });
    
    var  vegMean = ee.List([vegMean.get('b'),vegMean.get('g'),vegMean.get('r'),vegMean.get('nir'),vegMean.get('swir1'),vegMean.get('swir2')])
    
    var bareMean = lsuImageBS.reduceRegion({
      reducer: ee.Reducer.mean(),
      geometry: region,
      scale: resolution,
      tileScale: tileScale,
      maxPixels: 3e10
    });
    
    var  bareMean = ee.List([bareMean.get('b'),bareMean.get('g'),bareMean.get('r'),bareMean.get('nir'),bareMean.get('swir1'),bareMean.get('swir2')])
    
    // var lsuImage = lsuImage.unmix([waterMean, vegMean, bareMean]);
  
    var endmembers = ee.Feature(null,{"Water":waterMean})
    var endmembers = endmembers.set({"Veg":vegMean})
    var endmembers = endmembers.set({"Bare":bareMean})
    
    return endmembers
  }
  
  var endMembersFt = ee.Feature(genEndmembers(lsuImage));
  
  print('EndMembers: ',endmembers)
  
  var water = endMembersFt.get("Water")
  var veg = endMembersFt.get("Veg")
  var bare = endMembersFt.get("Bare")
  
  var waterMean = ee.List(water)
  var vegMean = ee.List(veg)
  var bareMean = ee.List(bare)
  
  var listVal = waterMean.get(0)
  // print(listVal)
  var waterMean = ee.Algorithms.If(listVal,waterMean,[0,0,0,0,0,0])
  // print("waterMean Check")
  // print(waterMean)
  
  var listVal = vegMean.get(0)
  // print(listVal)
  var vegMean = ee.Algorithms.If(listVal,vegMean,[0,0,0,0,0,0])
  // print("vegMean Check")
  // print(vegMean)
  
  var listVal = bareMean.get(0)
  // print(listVal)
  var bareMean = ee.Algorithms.If(listVal,bareMean,[0,0,0,0,0,0])
  // print("bareMean Check")
  // print(bareMean)
  
  var endmembers = ee.List([waterMean,vegMean,bareMean])
  
  // run the unmixing algorith
  var unmixedImage = lsuImage.unmix(endmembers,true,true);
  var bandNames2 = ['water', 'veg', 'bare'];
  var unmixedImage = unmixedImage.rename(bandNames2);
  var unmixedImage = unmixedImage
      .select(bandNames2);
      
  // print(unmixedImage);
  
  // Map.addLayer(unmixedImage,{min:0,max:100},"Unmixed")
  
  /////   STEP 3: IMAGE FRACTION CALCULATION
  
  // Sum to unity and calculate area fraction for each polygon.
  // This is completed in two steps to remove minus values 
  // thereby allowing fractions to sum to unity (b1>0<1).
  
  // Step 1: Assume first of all that all pixels contain minus values.
  // Reduce the image to get a one-band maximum value image.
  var maxValue = unmixedImage.reduce(ee.Reducer.max());
  //Map.addLayer(maxValue, {max: 13000}, 'Maximum value image');
  
  // Reduce the image to get a one-band minimum value image.
  var minValue = unmixedImage.reduce(ee.Reducer.min());
  //Map.addLayer(minValue, {max: 13000}, 'Minimum value image');
  
  // Subtract min from max as part of normalisation calculation:
  var rangeImage = maxValue.subtract(minValue);
  //Map.addLayer(rangeImage, {max: 13000}, 'Range'); 
  
  // Normalise the data to allow sum to unity per pixel normalisation = ((value - min)/(max - min)):
  var normImage = unmixedImage.subtract(minValue)
                  .divide(rangeImage);
  //Map.addLayer(normImage, {max: 13000}, 'Normalised'); 
  
  // Find the normalised sum of band fractions:
  var normSum = normImage.reduce(ee.Reducer.sum());
  //Map.addLayer(normSum, {max: 13000}, 'Normalised sum');
  
  // Find % fraction using normalised sum as total
  var fractionP = normImage.divide(normSum)
                  .multiply(100);
  //Map.addLayer(fractionP, {min:0, max: 100}, 'Fraction %'); 
  
  
  ///// STEP 4: BURN REGION DETECTION AND MASKING
  
  /*
    Confusion between burned pasture during the dry season and open water results
    in overestimation of the water fraction. To resolve this, a burn mask can be 
    applied to the study area as follows.
  */
  
  // Calculate the Burn Area Index (Chuvieco et al., 2002):
  
  var BAI = lsuImage.expression('float(1 / (((0.1-b("r"))*(0.1-b("r"))) + ((0.06-b("nir"))*(0.06-b("nir")))))');
  
  // Rename band:
  var BAIrenamed = BAI.select(
      ['constant'], // old names
      ['BAI',]      // new names
  );

  // Calculate metrics of central tendency for automatic thresholding (BAI):
  var BAIMean = BAIrenamed.reduceRegion({
    reducer: ee.Reducer.mean(),
    geometry: region,
    scale: resolution,
    maxPixels: 3e9,
    tileScale: tileScale
  }).values();
  
  // print('BAI Mean: ', BAIMean);
  var BAIstd = BAIrenamed.reduceRegion({
    reducer: ee.Reducer.stdDev(),
    geometry: region,
    scale: resolution,
    maxPixels: 3e9,
    tileScale: tileScale
  }).values();
  
  // print('BAI Std: ', BAIstd);
  
  // Calculate open water/bare sand fraction ratio where
  // 1 = pure open water and -1 = pure bare sand:
  
  // var OwBs = fractionP.expression('float((b("water")-b("bare"))/(b("water")+b("bare")))');
  var OwBs = fractionP.normalizedDifference(['water', 'bare'])
  // Map.addLayer(OwBs, {},  'OwBs_ratio');
  
  // Apply the threshold:
  // print("Threshold:",ee.Number(BAIMean.add(BAIstd)))
  
  var burnThresh4 = BAIrenamed.expression(
        "float(b('BAI') > (Mean1+Std1))",
        {
        'Mean1': ee.Number(BAIMean),
        'Std1': ee.Number(BAIstd)
  }).and(OwBs.gt(-0.5)).and(OwBs.lt(0.5)).and(fractionP.select('veg').lt(10));
  

  // Map.addLayer(burnThresh4,{},"BurnThresh") 
  
  
  ///// STEP 5: BURN MASK APPLICATION
  
  /*
    It should be noted that in the unmixing of multispectral imagery with band numbers of <= 8, soil and dry 
    bare sand are treated synonymously to prevent overfitting of the mixture model. Therefore, the following
    mention of bare sand associated with burned pasture is in relation to the exposure of bare sediments as 
    overlying vegtation is stripped away. In most instances, burned pasture appears high in the bare sand 
    fraction and low in water and veg fractions as underlying sediment is exposed. However, in freshly burned 
    regions that appear extremely dark, confusion occurs between water and bare sand. 
    This can be eliminated by subtractng the water fraction from itself and summing the subtracted value with 
    the bare sand fraction.
  */
  
  var burn_regions = burnThresh4;
  var burn_V = fractionP.updateMask(burn_regions);
  
  var waterUpdate = burn_V.select('water').subtract(burn_V.select('water'))
  var bareUpdate = burn_V.select('bare').add(burn_V.select('water'))
  
  var toMosaicWater = ee.ImageCollection([fractionP.select("water"), waterUpdate]);

  var mosaic2 = toMosaicWater.mosaic();
  var toMosaicBare = ee.ImageCollection([fractionP.select("bare"), bareUpdate]);
  var mosaic3 = toMosaicBare.mosaic(); 
  
  var burnRemoved = fractionP.addBands([mosaic2, mosaic3]);
  var burnRemovedFinal = burnRemoved.select(
      ['water_1', 'veg', 'bare_1'], // old names
      ['water', 'veg', 'bare',]               // new names
  );
  
  var fractionP = burnRemovedFinal;
  
  var fractionP = fractionP.toUint8()
  
  // Map.addLayer(fractionP, {min:0, max: 100}, 'Fraction % burn masked'); 
  
  // print(fractionP)
  
//   /////   STEP 6: THEMATIC CLASSIFICATION
  
  
// // Classification of image fractions where:
//         // 225 m2 = 25 % of pixel area
//         // 675 m2 = 75 % of pixel area
//         // 810 m2 = 90 % of pixel area
//         // 855 m2 = 95 % of pixel area
//         // 891 m2 = 99 % of pixel area
//         // 900 m2 = 100 % of pixel area 
// // Rule-based classification is as follows where:
//         // OW = Water 
//         // GV = Veg
//         // BS = Sand 
// // User can vary what is considered as the 'dominant' fraction by altering metres squared (m2).
// // Default is as follows:
//         // If Water and Veg cover >= 75% of a pixel 
//             // and Water alone is >= 25% but < 75%
//             // and Veg alone is >= 25% but < 75%
//             // then: Emergent Flooded Vegetation        (EFV)	(1)
//         // If Water and Sand cover >= 75% of a pixel 
//             // and Water alone is >= 25% but < 75%
//             // and Sand alone is >= 25% but < 75%
//             // then: Wet Bare Sand                      (WBS)	(2)
//         // If Veg and Sand cover >= 75% of a pixel 
//             // and Veg alone is >= 25% but < 75%
//             // and Sand alone is >= 25% but < 75%
//             // then: Dry Vegetation (sparse)            (DV)	(3)
//         // If Water alone covers >= 75% of a pixel 
//             // then: Open Water                         (OW1)	(4)
//         // If Veg alone covers >= 75% of a pixel 
//             // then: Dense Vegetation                   (GV1)	(5)
//         // If Sand alone covers >= 75% of a pixel 
//             // then: Dry Bare Sand / Barren             (BS1)	(6)	
//         // Else: The pixel is defined as Mixed/Other		(MIX)	(0)		

// Sum fraction images:

// OW summed with GV:
var WV = fractionP.select('water').add(fractionP.select('veg'));
//Map.addLayer(WV, {max: 13000}, 'water_veg');
// OW summed with BS:
var WB = fractionP.select('water').add(fractionP.select('bare'));
// Map.addLayer(WB, {max: 13000}, 'water_bare');
// GV summed with BS:
var VB = fractionP.select('veg').add(fractionP.select('bare'));
//Map.addLayer(VB, {max: 13000}, 'veg_bare');

// Thresholding to create conditional binary images:

//  set thresholds using expression function.
// Allows >= or <= to be set more easily and simplifies code: 
var WBS = fractionP.expression(
    //'(WB >= 75) & (water >= 25 & water < 75) & (bare >= 25 & bare < 75) & (slope <=2.5)', {
    '(WB >= 75) & (water >= 25 & water < 75) & (bare >= 25 & bare < 75) & (slope <=2.5)', { 
        'WB': WB,
        'water': fractionP.select('water'),
        'bare': fractionP.select('bare'),
        'HAND': hand.select('b1'),
        'slope': slope.select('slope')
});
var DV = fractionP.expression(
    '(((VB >= 75) & (veg >= 25 & veg < 75) & (bare >= 25 & bare < 75)) || ((WV >= 75) & (veg >= 25 & veg < 75) & (water >= 25 & water < 75) & (slope >2.5)))', {
        'VB': VB,
        'veg': fractionP.select('veg'),
        'bare': fractionP.select('bare'),
        'slope': slope.select('slope'),
        'WV': WV,
        'water': fractionP.select('water'),
        'HAND' : hand.select('b1')
});

var OW1 = fractionP.expression(
    '(water >= 75) & (slope <=2.5) & (HAND<=30)', {
        'water': fractionP.select('water'),
        'HAND': hand.select('b1'),
        'slope': slope.select('slope')
});


//// Mask out the DV and OW so we have only Dry and wet veg reamining; inundated veg is here \\\\
var EFV = fractionP.expression(
    //'(WV >= 75) & (veg >= 25 & veg < 75) & (water >= 25 & water < 75) & (slope <=0.5)', {
    '(WV >= 75) & (veg >= 25 & veg < 75) & (water >= 50 & water < 75) & (slope <=0.5)', {
        'WV': WV,
        'veg': fractionP.select('veg'),
        'water': fractionP.select('water'),
        'slope': slope.select('slope')

});

var GV1 = fractionP.expression(
    '(veg >= 75)', {
        'veg': fractionP.select('veg')
});


// NDWI mask based on 70th percentile
var GVNDII = ndiiComp.updateMask(GV1)
// print('GVNDII',GVNDII)
// Map.addLayer(GV1,{},"GV")
// Map.addLayer(GVNDII,{},"GVNDII")

////70th Percentile \\\\

var GVNDII_perc = GVNDII.reduceRegion({
    reducer: ee.Reducer.percentile([70]),
    geometry: geometry,
    scale: resolution,
    tileScale:tileScale,
    maxPixels:1e10
  });
// print(GVNDII_perc)
var GVNDII_perc = ee.Number(GVNDII_perc.get("nd"))
var AQ = GVNDII.gte(GVNDII_perc)
var AQ = ee.Image(ee.Algorithms.If(ee.Number(GVNDII_perc).gt(ee.Number(0.075)),AQ,ee.Image(0)))

// print("70th Percentile:",GVNDII_perc)
// Map.addLayer(AQ,{},"70th Percentile Threshold")

//// Error Logic Check \\\\
var GV2 = GV1.where(AQ.eq(1),0)

var AQ = AQ.unmask(0)

var AQ = AQ.updateMask(handSlopeMask)

// Map.addLayer(GV1,{},"GV1")

// Map.addLayer(GV2,{},"GV2")

var BS1 = fractionP.expression(
    '(bare >= 75)', {
        'bare': fractionP.select('bare')
});

  // Multiply conditional binary images by desired class value:
  var EFV_1 = ee.Image(EFV.multiply(1).rename("EFV"));
  var WBS_2 = ee.Image(WBS.multiply(2).rename("WBS"));
  var DV_3 = ee.Image(DV.multiply(3).rename("DV"));
  var OW1_4 = ee.Image(OW1.multiply(4).rename("OW"));
  var GV1_5 = ee.Image(GV2.multiply(5).rename("GV"));
  var BS1_6 = ee.Image(BS1.multiply(6).rename("BS"));
  var AQ_7 = ee.Image(AQ.multiply(7).rename("AQ"))
  
  // Map.addLayer(AQ_7,{},"AQ_7")
  
  var classImage2 = ee.Image(EFV_1.add(WBS_2).add(DV_3).add(OW1_4).add(GV1_5).add(BS1_6).add(AQ_7))

  // option to mask using information on forest cover, i.e.
  // TropWet unlikely to work over forestes areas
  
  // var dataset = ee.ImageCollection("COPERNICUS/Landcover/100m/Proba-V/Global")
  // var forestType = ee.Image(dataset.first().select('discrete_classification')).gte(111)
  // var forestType=forestType.updateMask(forestType)
  // print(forestType)
  
  // Map.addLayer(forestType.updateMask(forestType),{min:0,max:5},"forestType")

  // var forestMask = forestNonForest.neq(1)
  
  // var classImage2 = classImage2.where(forestType.eq(1),8)
  // Create a new multiband layer containing multiplied binary layers:
  // var imageUpdate3 = fractionP.addBands([EFV_1, WBS_2, DV_3, OW1_4, GV1_5, BS1_6]);
  //Map.addLayer(imageUpdate3, {max: 13000}, 'check');
  // Select and rename bands.
  // print(imageUpdate3)
  // var renamed = imageUpdate3.select(
  //     ['water_1', 'water_1_1', 'veg_1', 'water_2', 'veg_1_1', 'bare_1'], // old names
  //     ['EFV', 'WBS', 'DV', 'OW1', 'GV1', 'BS1',]               // new names
  // );
  
  // var renamed = imageUpdate3.select(['water_max','class'])
  
  // print(renamed)
  
  // // Sum all bands using image reducer to create a single band class image:
  // var classImage2 = renamed.reduce(ee.Reducer.sum());
  
  var vizParams3 = {
    min: 0,
    max: 8,
  };
  
  //Map.addLayer(classImage2, vizParams3, 'Class_Image');
  
  print(classImage2);
  
  // Ensure that NoData areas are exported as -9999 rather than 0
  
  var classImage3 = classImage2.unmask(0).clip(geometry);
  // print(classImage3)
  
  
  ///// STEP 5: DISPLAYING CLASS IMAGE AND LEGEND
  
  // Define a palette for the 18 distinct land cover classes.
  var classPalette = [
    'ffffff', // Cloud (Cloud, 0)
    '45AE90', // Emergent Flooded Veg (EFV, 1)
    '603114', // Wet Bare Sand (WBS, 2)
    'E5CB71', // Dry Sparse Vegetation (DV, 3)
    '2250CA', // Open Water (OW, 4)
    '328536', // Dense Vegetation (GV1, 5)
    'CBBF7C', // Bare Sand (BS1, 6)
    '00CCFF', // Aquatic Veg (AQ, 7)
    '000000', // Unsuitable Areas (Modis Mask, 8)
  ];
  
  // Specify the min and max labels and the color palette matching the labels.
  Map.addLayer(classImage3,
              {min: 0, max: 8.1, palette: classPalette},
              'Class_Image_Coloured', true);
  
  // set position of panel
  var legend = ui.Panel({
    style: {
      position: 'bottom-left',
      padding: '8px 15px'
    }
  });
   
  // Create legend title
  var legendTitle = ui.Label({
    value: 'Fraction Class',
    style: {
      fontWeight: 'bold',
      fontSize: '18px',
      margin: '0 0 4px 0',
      padding: '0'
      }
  });
   
  // Add the title to the panel
  legend.add(legendTitle);
   
  // Creates and styles 1 row of the legend.
  var makeRow = function(color, name) {
   
        // Create the label that is actually the colored box.
        var colorBox = ui.Label({
          style: {
            backgroundColor: '#' + color,
            // Use padding to give the box height and width.
            padding: '8px',
            margin: '0 0 4px 0' 
          }
        });
   
        // Create the label filled with the description text.
        var description = ui.Label({
          value: name,
          style: {margin: '0 0 4px 6px'}
        });
   
        // return the panel
        return ui.Panel({
          widgets: [colorBox, description],
          layout: ui.Panel.Layout.Flow('horizontal')
        });
  };
   
  //  Palette with the colors
  var palette =['ffffff', '45AE90','603114', 'E5CB71', '2250CA', '328536', 'CBBF7C','00CCFF','000000'];
   
  // name of the legend
  var names = ['Unclassified', 'Emergent Flooded Veg', 'Wet Bare Sediment', 'Dry Sparse Vegetation', 'Open Water', 'Dense Vegetation', 'Dry Bare Sediment','Aquatic Vegetation','Unsuitable'];
   
  // Add color and and names
  for (var i = 0; i < 9; i++) {
    legend.add(makeRow(palette[i], names[i]));
    }  
   
  // add legend to map (alternatively you can also print the legend to the console)
  Map.add(legend);
      
  
  ////
  //// EXPORTING OUTPUT PRODUCTS
  //
  
  // for export, convert to int, extract nir,swir,g
  var median_composite_scaled_export=median_composite_scaled.select('b','nir','swir1').toInt16()

  // Use to export fraction images in m2:
  Export.image.toDrive({
    image: fractionP,
    folder: '',
    description: 'Fractions_'+poly_ind,
    scale: resolution,
    region: region,
    maxPixels: 3e9
  }); 
  
  // export the thematic classification
  Export.image.toDrive({
    image: classImage3,
    folder: '',
    description: 'Classified_Output_'+poly_ind,
    scale: resolution,
    region: region,
    maxPixels: 3e9
  }); 
  
}