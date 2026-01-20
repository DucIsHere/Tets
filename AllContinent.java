package com.regenerationforrged.world.worldgen.cell.continent;

import com.regenerationforrged.world.worldgen.cell.Cell;
import com.regenerationforrged.world.worldgen.cell.CellPopulator;
import com.regenerationforrged.world.worldgen.cell.rivermap.Rivermap;

public interface Continent extends CellPopulator {
    float getEdgeValue(float x, float z);
    
    default float getLandValue(float x, float z) {
        return this.getEdgeValue(x, z);
    }
    
    long getNearestCenter(float x, float z);
    
    Rivermap getRivermap(int x, int z);
    
    default Rivermap getRivermap(Cell cell) {
        return this.getRivermap(cell.continentX, cell.continentZ);
    }
}

---

package com.regenerationforrged.world.worldgen.cell.continent;

import com.regenerationforrged.world.worldgen.cell.Cell;
import com.regenerationforrged.world.worldgen.cell.CellPopulator;
import com.regenerationforrged.world.worldgen.noise.NoiseUtil;
import com.regenerationforrged.world.worldgen.noise.function.Interpolation;

public class ContinentLerper2 implements CellPopulator {
	private CellPopulator lower;
    private CellPopulator upper;
    private Interpolation interpolation;
    private float blendLower;
    private float blendUpper;
    private float blendRange;
    
    public ContinentLerper2(CellPopulator lower, CellPopulator upper, float min, float max) {
        this(lower, upper, min, max, Interpolation.LINEAR);
    }
    
    public ContinentLerper2(CellPopulator lower, CellPopulator upper, float min, float max, Interpolation interpolation) {
        this.lower = lower;
        this.upper = upper;
        this.interpolation = interpolation;
        this.blendLower = min;
        this.blendUpper = max;
        this.blendRange = this.blendUpper - this.blendLower;
    }
    
    @Override
    public void apply(Cell cell, float x, float y) {
        if (cell.continentEdge < this.blendLower) {
            this.lower.apply(cell, x, y);
            return;
        }
        if (cell.continentEdge > this.blendUpper) {
            this.upper.apply(cell, x, y);
            return;
        }
        float alpha = this.interpolation.apply((cell.continentEdge - this.blendLower) / this.blendRange);
        this.lower.apply(cell, x, y);
        float lowerVal = cell.height;
        this.upper.apply(cell, x, y);
        float upperVal = cell.height;
        cell.height = NoiseUtil.lerp(lowerVal, upperVal, alpha);
    }
}

---

package com.regenerationforrged.world.worldgen.cell.continent;

import com.regenerationforrged.world.worldgen.cell.Cell;
import com.regenerationforrged.world.worldgen.cell.CellPopulator;
import com.regenerationforrged.world.worldgen.noise.NoiseUtil;
import com.regenerationforrged.world.worldgen.noise.function.Interpolation;

public class ContinentLerper3 implements CellPopulator {
    private CellPopulator lower;
    private CellPopulator middle;
    private CellPopulator upper;
    private Interpolation interpolation;
    private float midpoint;
    private float blendLower;
    private float blendUpper;
    private float lowerRange;
    private float upperRange;
    
    public ContinentLerper3(CellPopulator lower, CellPopulator middle, CellPopulator upper, float min, float mid, float max) {
        this(lower, middle, upper, min, mid, max, Interpolation.CURVE3);
    }
    
    public ContinentLerper3(CellPopulator lower, CellPopulator middle, CellPopulator upper, float min, float mid, float max, Interpolation interpolation) {
        this.lower = lower;
        this.upper = upper;
        this.middle = middle;
        this.interpolation = interpolation;
        this.midpoint = mid;
        this.blendLower = min;
        this.blendUpper = max;
        this.lowerRange = this.midpoint - this.blendLower;
        this.upperRange = this.blendUpper - this.midpoint;
    }
    
    @Override
    public void apply(Cell cell, float x, float y) {
        float select = cell.continentEdge;
        if (select < this.blendLower) {
            this.lower.apply(cell, x, y);
            return;
        }
        if (select > this.blendUpper) {
            this.upper.apply(cell, x, y);
            return;
        }
        if (select < this.midpoint) {
            float alpha = this.interpolation.apply((select - this.blendLower) / this.lowerRange);
            this.lower.apply(cell, x, y);
            float lowerVal = cell.height;
            this.middle.apply(cell, x, y);
            cell.height = NoiseUtil.lerp(lowerVal, cell.height, alpha);
        } else {
            float alpha = this.interpolation.apply((select - this.midpoint) / this.upperRange);
            this.middle.apply(cell, x, y);
            float lowerVal = cell.height;
            this.upper.apply(cell, x, y);
            cell.height = NoiseUtil.lerp(lowerVal, cell.height, alpha);
        }
    }
}

---

package com.regenerationforrged.world.worldgen;

import org.jetbrains.annotation.Nullable;

import net.minecraft.core.HolderGetter;

public class GeneratorContext {
  public Seed seed;
  public Levels levels;
  public Preset preset;
  public HolderGetter<Noise> noiseLookup;
  public TitleGenerator generator;
  @Nullable
  public TitleCache cache;
  public WorldLookup Lookup;

  public GeneratorContext(Preset preset, HolderGetter<Noise> noiseLookup, int seed, int titleSize, int titleBorder, int batchCount, @Nullable TitleCache cache) {
    this.preset = preset;
    this.noiseLookup = noiseLookup;
    this.seed = new Seed(seed);
    this.levels = new Levels(preset.world().properties.terrianScaler(), preset.world().properties.seaLevel);
    this.cache = cache;
    
  }
}

---

package com.regenerationforrged.world.worldgen.cell.continent;

import com.regenerationforrged.world.worldgen.cell.Cell;
import com.regenerationforrged.world.worldgen.cell.CellPopulator;
import com.regenerationforrged.world.worldgen.noise.NoiseUtil;

public record MushroomIslandPopulator(CellPopulator ocean, CellPopulator coast, CellPopulator land, float coastPoint, float inlandPoint) implements CellPopulator {
	public static final float DEFAULT_INLAND_POINT = 0.0F;
	public static final float DEFAULT_COAST_POINT = DEFAULT_INLAND_POINT + 0.1F;
	
	@Override
	public void apply(Cell cell, float x, float z) {
		if(cell.continentEdge < this.inlandPoint) {
			this.land.apply(cell, x, z);
			return;
		}
		
		if(cell.continentEdge < this.coastPoint) {
			this.coast.apply(cell, x, z);
			float coastHeight = cell.height;
			
			this.land.apply(cell, x, z);
			float landHeight = cell.height;
			
			cell.height = NoiseUtil.lerp(landHeight, coastHeight, NoiseUtil.map(cell.continentEdge, 0.0F, 1.0F, this.inlandPoint, this.coastPoint));
			return;
		}
		
		this.ocean.apply(cell, x, z);
	}
}

---

package com.regenerationforrged.world.worldgen.cell.continent;

public interface SimpleContinent extends Continent {
  floatEdgeValue(float x, float z);

  default float getDistanceToEdge(int cx, int cz, float dx, float dz) {
    return 1.0F;
  }

  default float getDistanceToOcean(int cx, int cz, float dx, float dy) {
    return 1.0F;
  }
}

---

package com.regenerationforrged.world.worldgen.cell.advaced;

import com.regenerationforrged.data.worldgen.preset.settings.WorldSettings;
import com.regenerationforrged.world.worldgen.GeneratorContext;
import com.regenerationforrged.world.worldgen.cell.continent.SimpleContinent;
import com.regenerationforrged.world.worldgen.cell.continent.simple.SimpleRiverGenerator;
import com.regenerationforrged.world.worldgen.cell.heightmap.ControlPoints;
import com.regenerationforrged.world.worldgen.cell.rivermap.RiverCache;
import com.regenerationforrged.world.worldgen.noise.NoiseUtil;
import com.regenerationforrged.world.worldgen.util.PosUtil;
import com.regenerationforrged.world.worldgen.util.Seed;

public abstract class AbstractContinent implements SimpleContinent {
    protected int seed;
    protected int skippingSeed;
    protected int continentScale;
    protected float jitter;
    protected boolean hasSkipping;
    protected float skipThreshold;
    protected RiverCache riverCache;
    protected ControlPoints controlPoints;
    
    public AbstractContinent(Seed seed, GeneratorContext context) {
        WorldSettings settings = context.preset.world();
        this.seed = seed.next();
        this.skippingSeed = seed.next();
        this.continentScale = settings.continent.continentScale;
        this.jitter = settings.continent.continentJitter;
        this.skipThreshold = settings.continent.continentSkipping;
        this.hasSkipping = (this.skipThreshold > 0.0F);
        this.controlPoints = ControlPoints.make(settings.controlPoints);
        this.riverCache = new RiverCache(new SimpleRiverGenerator(this, context));
    }
    
    @Override
    public float getDistanceToOcean(int cx, int cz, float dx, float dz) {
        float high = this.getDistanceToEdge(cx, cz, dx, dz);
        float low = 0.0F;
        for (int i = 0; i < 50; ++i) {
            float mid = (low + high) / 2.0F;
            float x = cx + dx * mid;
            float z = cz + dz * mid;
            float edge = this.getEdgeValue(x, z);
            if (edge > this.controlPoints.shallowOcean()) {
                low = mid;
            } else {
                high = mid;
            }
            if (high - low < 10.0F) {
                break;
            }
        }
        return high;
    }
    
    @Override
    public float getDistanceToEdge(int cx, int cz, float dx, float dz) {
        float distance = (float)(this.continentScale * 4);
        for (int i = 0; i < 10; ++i) {
            float x = cx + dx * distance;
            float z = cz + dz * distance;
            long centerPos = this.getNearestCenter(x, z);
            int conX = PosUtil.unpackLeft(centerPos);
            int conZ = PosUtil.unpackRight(centerPos);
            distance += distance;
            if (conX != cx || conZ != cz) {
                float low = 0.0f;
                float high = distance;
                for (int j = 0; j < 50; ++j) {
                    float mid = (low + high) / 2.0F;
                    float px = cx + dx * mid;
                    float pz = cz + dz * mid;
                    centerPos = this.getNearestCenter(px, pz);
                    conX = PosUtil.unpackLeft(centerPos);
                    conZ = PosUtil.unpackRight(centerPos);
                    if (conX == cx && conZ == cz) {
                        low = mid;
                    } else {
                        high = mid;
                    }
                    if (high - low < 50.0F) {
                        break;
                    }
                }
                return high;
            }
        }
        return distance;
    }

    protected boolean isDefaultContinent(int cellX, int cellY) {
        return cellX == 0 && cellY == 0;
    }
    
    protected boolean shouldSkip(int cellX, int cellY) {
        if (this.hasSkipping && !this.isDefaultContinent(cellX, cellY)) {
            float skipValue = getCellValue(this.skippingSeed, cellX, cellY);
            return skipValue < this.skipThreshold;
        }
        return false;
    }
    
    protected static float getCellValue(int seed, int cellX, int cellY) {
        return 0.5F + NoiseUtil.valCoord2D(seed, cellX, cellY) * 0.5F;
    }
}

---


package com.regenerationforrged.world.worldgen.cell.continent.advanced;

import com.regenerationforrged.concurrent.Resource;
import com.regenerationforrged.data.worldgen.preset.settings.WorldSettings;
import com.regenerationforrged.world.worldgen.GeneratorContext;
import com.regenerationforrged.world.worldgen.cell.Cell;
import com.regenerationforrged.world.worldgen.cell.continent.SimpleContinent;
import com.regenerationforrged.world.worldgen.cell.rivermap.Rivermap;
import com.regenerationforrged.world.worldgen.noise.NoiseUtil;
import com.regenerationforrged.world.worldgen.noise.NoiseUtil.Vec2f;
import com.regenerationforrged.world.worldgen.noise.domain.Domain;
import com.regenerationforrged.world.worldgen.noise.domain.Domains;
import com.regenerationforrged.world.worldgen.noise.module.Line;
import com.regenerationforrged.world.worldgen.noise.module.Noise;
import com.regenerationforrged.world.worldgen.noise.module.Noises;
import com.regenerationforrged.world.worldgen.util.PosUtil;
import com.regenerationforrged.world.worldgen.util.Seed;

public class AdvancedContinentGenerator extends AbstractContinent implements SimpleContinent {
    protected static float CENTER_CORRECTION = 0.35F;
    protected float frequency;
    protected float variance;
    protected int varianceSeed;
    protected Domain warp;
    protected Noise cliffNoise;
    protected Noise bayNoise;
    
    public AdvancedContinentGenerator(Seed seed, GeneratorContext context) {
        super(seed, context);
        WorldSettings settings = context.preset.world();
        int tectonicScale = settings.continent.continentScale * 4;
        this.frequency = 1.0F / tectonicScale;
        this.varianceSeed = seed.next();
        this.variance = settings.continent.continentSizeVariance;
        this.warp = this.createWarp(seed, tectonicScale, settings.continent);
        
        float frequency = 1.0F / this.frequency;
        
        Noise cliffNoise = Noises.simplex2(seed.next(), this.continentScale / 2, 2);
        cliffNoise = Noises.clamp(cliffNoise, 0.1F, 0.25F);
        cliffNoise = Noises.map(cliffNoise, 0.0F, 1.0F);
        cliffNoise = Noises.frequency(cliffNoise, frequency);
        this.cliffNoise = cliffNoise;

        Noise bayNoise = Noises.simplex(seed.next(), 100, 1);
        bayNoise = Noises.mul(bayNoise, 0.1F);
        bayNoise = Noises.add(bayNoise, 0.9F);
        bayNoise = Noises.frequency(bayNoise, frequency);
        this.bayNoise = bayNoise;
    }
    
    @Override
    public void apply(Cell cell, float x, float y) {
        float wx = this.warp.getX(x, y, 0);
        float wy = this.warp.getZ(x, y, 0);
        x = wx * this.frequency;
        y = wy * this.frequency;
        int xi = NoiseUtil.floor(x);
        int yi = NoiseUtil.floor(y);
        int cellX = xi;
        int cellY = yi;
        float cellPointX = x;
        float cellPointY = y;
        float nearest = Float.MAX_VALUE;
        for (int cy = yi - 1; cy <= yi + 1; ++cy) {
            for (int cx = xi - 1; cx <= xi + 1; ++cx) {
                Vec2f vec = NoiseUtil.cell(this.seed, cx, cy);
                float px = cx + vec.x() * this.jitter;
                float py = cy + vec.y() * this.jitter;
                float dist2 = Line.distSq(x, y, px, py);
                if (dist2 < nearest) {
                    cellPointX = px;
                    cellPointY = py;
                    cellX = cx;
                    cellY = cy;
                    nearest = dist2;
                }
            }
        }
        nearest = Float.MAX_VALUE;
        float sumX = 0.0f;
        float sumY = 0.0f;
        for (int cy2 = cellY - 1; cy2 <= cellY + 1; ++cy2) {
            for (int cx2 = cellX - 1; cx2 <= cellX + 1; ++cx2) {
                if (cx2 != cellX || cy2 != cellY) {
                    Vec2f vec2 = NoiseUtil.cell(this.seed, cx2, cy2);
                    float px2 = cx2 + vec2.x() * this.jitter;
                    float py2 = cy2 + vec2.y() * this.jitter;
                    float dist3 = getDistance(x, y, cellPointX, cellPointY, px2, py2);
                    sumX += px2;
                    sumY += py2;
                    if (dist3 < nearest) {
                        nearest = dist3;
                    }
                }
            }
        }
        if (this.shouldSkip(cellX, cellY)) {
            return;
        }
        cell.continentId = AbstractContinent.getCellValue(this.seed, cellX, cellY);
        cell.continentEdge = this.getDistanceValue(x, y, cellX, cellY, nearest);
        cell.continentX = this.getCorrectedContinentCenter(cellPointX, sumX / 8.0F);
        cell.continentZ = this.getCorrectedContinentCenter(cellPointY, sumY / 8.0F);
    }
    
    @Override
    public float getEdgeValue(float x, float z) {
        try (Resource<Cell> resource = Cell.getResource()) {
            Cell cell = resource.get();
            this.apply(cell, x, z);
            return cell.continentEdge;
        }
    }
    
    @Override
    public long getNearestCenter(float x, float z) {
        try (Resource<Cell> resource = Cell.getResource()) {
            Cell cell = resource.get();
            this.apply(cell, x, z);
            return PosUtil.pack(cell.continentX, cell.continentZ);
        }
    }
    
    @Override
    public Rivermap getRivermap(int x, int z) {
        return this.riverCache.getRivers(x, z);
    }
    
    protected Domain createWarp(Seed seed, int tectonicScale, WorldSettings.Continent continent) {
        int warpScale = NoiseUtil.round(tectonicScale * 0.225F);
        float strength = NoiseUtil.round(tectonicScale * 0.33F);
        return Domains.domain(
        	Noises.perlin2(seed.next(), warpScale, continent.continentNoiseOctaves, continent.continentNoiseLacunarity, continent.continentNoiseGain), 
        	Noises.perlin2(seed.next(), warpScale, continent.continentNoiseOctaves, continent.continentNoiseLacunarity, continent.continentNoiseGain), 
        	Noises.constant(strength)
        );
    }
    
    protected float getDistanceValue(float x, float y, int cellX, int cellY, float distance) {
        distance = this.getVariedDistanceValue(cellX, cellY, distance);
        distance = NoiseUtil.sqrt(distance);
        distance = NoiseUtil.map(distance, 0.05f, 0.25f, 0.2f);
        distance = this.getCoastalDistanceValue(x, y, distance);
        if (distance < this.controlPoints.inland() && distance >= this.controlPoints.shallowOcean()) {
            distance = this.getCoastalDistanceValue(x, y, distance);
        }
        return distance;
    }
    
    protected float getVariedDistanceValue(int cellX, int cellY, float distance) {
        if (this.variance > 0.0f && !this.isDefaultContinent(cellX, cellY)) {
            float sizeValue = AbstractContinent.getCellValue(this.varianceSeed, cellX, cellY);
            float sizeModifier = NoiseUtil.map(sizeValue, 0.0f, this.variance, this.variance);
            distance *= sizeModifier;
        }
        return distance;
    }
    
    protected float getCoastalDistanceValue(float x, float y, float distance) {
        if (distance > this.controlPoints.shallowOcean() && distance < this.controlPoints.inland()) {
            float alpha = distance / this.controlPoints.inland();
            float cliff = this.cliffNoise.compute(x, y, 0);
            distance = NoiseUtil.lerp(distance * cliff, distance, alpha);
            if (distance < this.controlPoints.shallowOcean()) {
                distance = this.controlPoints.shallowOcean() * this.bayNoise.compute(x, y, 0);
            }
        }
        return distance;
    }
    
    protected int getCorrectedContinentCenter(float point, float average) {
        point = NoiseUtil.lerp(point, average, 0.35f) / this.frequency;
        return (int)point;
    }
    
    protected static float midPoint(float a, float b) {
        return (a + b) * 0.5F;
    }
    
    protected static float getDistance(float x, float y, float ax, float ay, float bx, float by) {
        float mx = midPoint(ax, bx);
        float my = midPoint(ay, by);
        float dx = bx - ax;
        float dy = by - ay;
        float nx = -dy;
        float ny = dx;
        return getDistance2Line(x, y, mx, my, mx + nx, my + ny);
    }
    
    protected static float getDistance2Line(float x, float y, float ax, float ay, float bx, float by) {
        float dx = bx - ax;
        float dy = by - ay;
        float v = (x - ax) * dx + (y - ay) * dy;
        v /= dx * dx + dy * dy;
        float ox = ax + dx * v;
        float oy = ay + dy * v;
        return Line.distSq(x, y, ox, oy);
    }
}

---

package com.regenerationforrged.world.worldgen.cell.continent.fancy;

import java.util.Random;

import com.regenerationforrged.world.worldgen.GeneratorContext;
import com.regenerationforrged.world.worldgen.cell.heightmap.ControlPoints;
import com.regenerationforrged.world.worldgen.cell.rivermap.RiverGenerator;
import com.regenerationforrged.world.worldgen.cell.rivermap.Rivermap;
import com.regenerationforrged.world.worldgen.noise.NoiseUtil;
import com.regenerationforrged.world.worldgen.noise.NoiseUtil.Vec2f;
import com.regenerationforrged.world.worldgen.util.PosUtil;

public class FancyContinent implements RiverGenerator {
	private Island[] islands;
	private FancyRiverGenerator riverGenerator;

	public FancyContinent(int seed, int nodes, float radius, GeneratorContext context, FancyContinentGenerator continent) {
		ControlPoints controlPoints = ControlPoints.make(context.preset.world().controlPoints);
		this.islands = generateIslands(controlPoints, 3, nodes, radius, new Random(seed));
		this.riverGenerator = new FancyRiverGenerator(continent, context);
	}

	public float getEdgeValue(float x, float y, int seed) {
		float value = 0.0F;
		for (Island island : this.islands) {
			float v = island.getEdgeValue(x, y);
			value = Math.max(v, value);
		}
		return process(value);
	}

	public Island getMain() {
		return this.islands[0];
	}

	public Island[] getIslands() {
		return this.islands;
	}

	public long getMin() {
		float x = Float.MAX_VALUE;
		float z = Float.MAX_VALUE;
		for (Island island : this.islands) {
			x = Math.min(x, island.getMin().x());
			z = Math.min(z, island.getMin().y());
		}
		return PosUtil.packf(x, z);
	}

	public long getMax() {
		float x = Float.MIN_VALUE;
		float z = Float.MIN_VALUE;
		for (Island island : this.islands) {
			x = Math.max(x, island.getMin().x());
			z = Math.max(z, island.getMin().y());
		}
		return PosUtil.packf(x, z);
	}

	public float getLandValue(float x, float y) {
		float value = 0.0F;
		for (Island island : this.islands) {
			float v = island.getLandValue(x, y);
			value = Math.max(v, value);
		}
		return value;
	}

	public long getValueId(float x, float y) {
		int id = -1;
		float value = 0.0F;
		for (Island island : this.islands) {
			float v = island.getEdgeValue(x, y);
			if (v > value) {
				value = v;
				id = island.getId();
			}
			value = Math.max(v, value);
		}
		return PosUtil.packMix(id, value);
	}

	@Override
	public Rivermap generateRivers(int x, int z, long id) {
		return this.riverGenerator.generateRivers(x, z, id);
	}

	private static float process(float value) {
		return value;
	}

	private static Island[] generateIslands(ControlPoints controlPoints, int islandCount, int nodeCount, float radius, Random random) {
		int dirs = 4;
		Island main = generate(0, controlPoints, nodeCount, radius, random);
		Island[] islands = new Island[1 + islandCount * dirs];
		islands[0] = main;
		int i = 1;
		float yawStep = 1.0F / dirs * 6.2831855F;
		for (int dir = 0; dir < dirs; ++dir) {
			Island previous = main;
			int nCount = Math.max(2, nodeCount - 1);
			float r = radius * 0.5F;
			float yaw = yawStep * dir + random.nextFloat() * yawStep;
			for (int island = 0; island < islandCount; ++island) {
				Island next = generate(i, controlPoints, nCount, r, random);
				float y = yaw + nextFloat(random, -0.2F, 0.2F);
				float distance = previous.radius();
				float dx = NoiseUtil.sin(y * 6.2831855F) * distance;
				float dz = NoiseUtil.cos(y * 6.2831855F) * distance;
				float ox = previous.getCenter().x() + dx;
				float oy = previous.getCenter().y() + dz;
				next.translate(new Vec2f(ox, oy));
				nCount = Math.max(2, nCount - 1);
				r *= 0.8F;
				islands[i++] = next;
				previous = next;
			}
		}
		return islands;
	}

	private static Island generate(int id, ControlPoints controlPoints, int nodes, float radius, Random random) {
		float minScale = 0.75F;
		float maxScale = 2.5F;
		float minLen = radius * 1.5F;
		float maxLen = radius * 3.5F;
		float maxYaw = 1.5707964F;
		float minYaw = -maxYaw;
		Segment[] segments = new Segment[nodes - 1];
		Vec2f pointA = new Vec2f(0.0F, 0.0F);
		float scaleA = nextFloat(random, minScale, maxScale);
		float previousYaw = nextFloat(random, 0.0F, 6.2831855F);
		for (int i = 0; i < segments.length; ++i) {
			float length = nextFloat(random, minLen, maxLen);
			float yaw = previousYaw + nextFloat(random, minYaw, maxYaw);
			float dx = NoiseUtil.sin(yaw) * length;
			float dz = NoiseUtil.cos(yaw) * length;
			Vec2f pointB = new Vec2f(pointA.x() + dx, pointA.y() + dz);
			float scaleB = nextFloat(random, minScale, maxScale);
			segments[i] = new Segment(pointA, pointB, scaleA, scaleB);
			previousYaw = yaw;
			pointA = pointB;
			scaleA = scaleB;
		}
		return new Island(id, segments, controlPoints, radius * 3.0f, radius * 1.25F, radius, radius * 0.975F);
	}
	
	public static float nextFloat(Random random, float min, float max) {
		return min + random.nextFloat() * (max - min);
	}
}

---

package com.regenerationforrged.world.worldgen.cell.continent.fancy;

import com.regenerationforrged.data.worldgen.preset.settings.WorldSettings;
import com.regenerationforrged.world.worldgen.GeneratorContext;
import com.regenerationforrged.world.worldgen.cell.Cell;
import com.regenerationforrged.world.worldgen.cell.continent.Continent;
import com.regenerationforrged.world.worldgen.cell.rivermap.RiverCache;
import com.regenerationforrged.world.worldgen.cell.rivermap.Rivermap;
import com.regenerationforrged.world.worldgen.noise.NoiseUtil;
import com.regenerationforrged.world.worldgen.noise.domain.Domain;
import com.regenerationforrged.world.worldgen.noise.domain.Domains;
import com.regenerationforrged.world.worldgen.util.PosUtil;
import com.regenerationforrged.world.worldgen.util.Seed;

public class FancyContinentGenerator implements Continent {
    private float frequency;
    private Domain warp;
    private FancyContinent source;
    private RiverCache riverCache;
    
    public FancyContinentGenerator(Seed seed, GeneratorContext context) {
        WorldSettings settings = context.preset.world();
        int warpScale = settings.continent.continentScale / 2;
        float warpStrength = warpScale * 0.4F;
        this.source = new FancyContinent(seed.next(), 4, 0.2F, context, this);
        this.frequency = 1.0F / settings.continent.continentScale;
        this.riverCache = new RiverCache(this.source);
        
        Domain warp = Domains.domainSimplex(seed.next(), warpScale, 2, warpStrength);
        warp = Domains.add(warp, Domains.domainPerlin(seed.next(), 80, 2, 40.0F)); 
        warp = Domains.add(warp, Domains.domainPerlin(seed.next(), 20, 1, 15.0F));
        this.warp = warp;
    }
    
    public FancyContinent getSource() {
        return this.source;
    }
    
    @Override
    public Rivermap getRivermap(int x, int y) {
        return this.riverCache.getRivers(x, y);
    }
    
    @Override
    public float getEdgeValue(float x, float y) {
        float px = this.warp.getX(x, y, 0);
        float py = this.warp.getZ(x, y, 0);
        px *= this.frequency;
        py *= this.frequency;
        return this.source.getEdgeValue(px, py, 0);
    }
    
    @Override
    public float getLandValue(float x, float y) {
        float px = this.warp.getX(x, y, 0);
        float py = this.warp.getZ(x, y, 0);
        px *= this.frequency;
        py *= this.frequency;
        float value = this.source.getLandValue(px, py);
        return NoiseUtil.map(value, 0.2F, 0.4F, 0.2F);
    }
    
    @Override
    public long getNearestCenter(float x, float y) {
        long min = this.source.getMin();
        long max = this.source.getMax();
        float width = PosUtil.unpackLeftf(max) - PosUtil.unpackLeftf(min);
        float height = PosUtil.unpackRightf(max) - PosUtil.unpackRightf(min);
        float cx = width * 0.5F;
        float cz = height * 0.5F;
        int centerX = (int)(cx / this.frequency);
        int centerZ = (int)(cz / this.frequency);
        return PosUtil.pack(centerX, centerZ);
    }
    
    @Override
    public void apply(Cell cell, float x, float y) {
        cell.continentX = 0;
        cell.continentZ = 0;
        cell.continentId = 0.0F;
        cell.continentEdge = this.getEdgeValue(x, y);
    }
}

---

package com.regenerationforrged.world.worldgen.cell.continent.fancy;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import com.regenerationforrged.world.worldgen.GeneratorContext;
import com.regenerationforrged.world.worldgen.cell.rivermap.Rivermap;
import com.regenerationforrged.world.worldgen.cell.rivermap.gen.GenWarp;
import com.regenerationforrged.world.worldgen.cell.rivermap.river.BaseRiverGenerator;
import com.regenerationforrged.world.worldgen.cell.rivermap.river.Network;
import com.regenerationforrged.world.worldgen.cell.rivermap.river.River;
import com.regenerationforrged.world.worldgen.cell.rivermap.river.RiverCarver;
import com.regenerationforrged.world.worldgen.cell.rivermap.river.RiverConfig;
import com.regenerationforrged.world.worldgen.cell.rivermap.river.RiverWarp;
import com.regenerationforrged.world.worldgen.noise.NoiseUtil;
import com.regenerationforrged.world.worldgen.noise.NoiseUtil.Vec2f;
import com.regenerationforrged.world.worldgen.util.PosUtil;
import com.regenerationforrged.world.worldgen.util.Variance;

public class FancyRiverGenerator extends BaseRiverGenerator<FancyContinentGenerator> {
	private static final Variance MAIN_PADDING = Variance.of(0.05F, 0.1F);
	private static final Variance MAIN_JITTER = Variance.of(-0.2F, 0.4F);
	private float freq;

	public FancyRiverGenerator(FancyContinentGenerator continent, GeneratorContext context) {
		super(continent, context);
		this.freq = 1.0F / context.preset.world().continent.continentScale;
	}

	@Override
	public Rivermap generateRivers(int x, int z, long id) {
		Random random = new Random(id + this.seed);
		GenWarp warp = GenWarp.make((int) id, this.continentScale);
		List<Network> networks = new ArrayList<>(32);
		List<Network.Builder> roots = new ArrayList<>(16);
		for (Island island : ((FancyContinentGenerator) this.continent).getSource().getIslands()) {
			this.generateRoots(((FancyContinentGenerator) this.continent).getSource(), island, random, warp, roots);
			for (Network.Builder river : roots) {
				networks.add(river.build());
			}
			roots.clear();
		}
		return new Rivermap(x, z, networks.toArray(Network[]::new), warp);
	}

	private void generateRoots(FancyContinent continent, Island island, Random random, GenWarp warp, List<Network.Builder> roots) {
		Segment[] segments = island.getSegments();
		int lineCount = Math.max(1, 8 - island.getId());
		int endCount = Math.max(4, 12 - island.getId());
		for (int i = 0; i < segments.length; ++i) {
			boolean end = i == 0 || i == segments.length - 1;
			Segment segment = segments[i];
			int riverCount = end ? (lineCount - 1) : lineCount;
			this.collectSegmentRoots(continent, island, segment, riverCount, random, warp, roots);
		}
		Segment first = segments[0];
		this.collectPointRoots(continent, island, first.a, first.scaleA, endCount, random, warp, roots);
		Segment last = segments[segments.length - 1];
		this.collectPointRoots(continent, island, last.b, last.scaleB, endCount, random, warp, roots);
	}

	private void collectSegmentRoots(FancyContinent continent, Island island, Segment segment, int count, Random random, GenWarp warp, List<Network.Builder> roots) {
		float dx = segment.dx;
		float dy = segment.dy;
		float nx = dy / segment.length;
		float ny = -dx / segment.length;
		float stepSize = 1.0F / (count + 2);
		for (int i = 0; i < count; ++i) {
			float progress = stepSize + stepSize * i;
			if (progress > 1.0F) {
				return;
			}
			float startX = segment.a.x() + dx * progress;
			float startZ = segment.a.y() + dy * progress;
			float radiusScale = NoiseUtil	.lerp(segment.scaleA, segment.scaleB, progress);
			float radius = island.coast() * radiusScale;
			int dir = random.nextBoolean() ? -1 : 1;
			float dirX = nx * dir + FancyRiverGenerator.MAIN_JITTER.next(random);
			float dirZ = ny * dir + FancyRiverGenerator.MAIN_JITTER.next(random);
			float scale = getExtendScale(island.getId(), startX, startZ, dirX, dirZ, radius, continent);
			if (scale != 0.0f) {
				float startPad = FancyRiverGenerator.MAIN_PADDING.next(random);
				float x1 = startX + dir * dirX * radius * startPad;
				float y1 = startZ + dir * dirZ * radius * startPad;
				float x2 = startX + dirX * radius * scale;
				float y2 = startZ + dirZ * radius * scale;
				this.addRoot(x1, y1, x2, y2, this.main, random, warp, roots);
			}
		}
	}

	private void collectPointRoots(FancyContinent continent, Island island, Vec2f vec, float radiusScale, int count, Random random, GenWarp warp, List<Network.Builder> roots) {
		float yawStep = 6.2831855f / count;
		float radius = island.coast() * radiusScale;
		for (int i = 0; i < count; ++i) {
			float yaw = yawStep * i;
			float dx = NoiseUtil.cos(yaw);
			float dz = NoiseUtil.sin(yaw);
			float scale = getExtendScale(island.getId(), vec.x(), vec.y(), dx, dz, radius, continent);
			if (scale != 0.0F) {
				float startPad = FancyRiverGenerator.MAIN_PADDING.next(random);
				float startX = vec.x() + dx * startPad * radius;
				float startZ = vec.y() + dz * startPad * radius;
				float endX = vec.x() + dx * radius * scale;
				float endZ = vec.y() + dz * radius * scale;
				if (continent.getEdgeValue(endX, endZ, seed) <= 0.1F) {
					this.addRoot(startX, startZ, endX, endZ, this.main, random, warp, roots);
				}
			}
		}
	}

	private void addRoot(float x1, float z1, float x2, float z2, RiverConfig config, Random random, GenWarp warp, List<Network.Builder> roots) {
		River river = new River(x1 / this.freq, z1 / this.freq, x2 / this.freq, z2 / this.freq);
		if (this.riverOverlaps(river, null, roots)) {
			return;
		}
		RiverCarver.Settings settings = BaseRiverGenerator.creatSettings(random);
		settings.fadeIn = config.fade;
		settings.valleySize = 275.0F * River.FORK_VALLEY.next(random);
		RiverWarp riverWarp = RiverWarp.create(0.1F, 0.85F, random);
		RiverCarver carver = new RiverCarver(river, riverWarp, config, settings, this.levels);
		Network.Builder network = Network.builder(carver);
		roots.add(network);
		this.generateForks(network, River.FORK_SPACING, this.fork, random, warp, roots, 0);
		this.generateWetlands(network, random);
	}

	private static float getExtendScale(int islandId, float startX, float startZ, float dx, float dz, float radius, FancyContinent continent) {
		float scale = 1.0F;
		for (int i = 0; i < 25; ++i) {
			float x = startX + dx * radius * scale;
			float z = startZ + dz * radius * scale;
			long packed = continent.getValueId(x, z);
			if (PosUtil.unpackLeft(packed) != islandId) {
				return 0.0F;
			}
			if (PosUtil.unpackRightf(packed) < 0.1F) {
				return scale;
			}
			scale += 0.075F;
		}
		return 0.0F;
	}
}

---

package com.regenerationforrged.world.worldgen.cell.continent.fancy;

import com.regenerationforrged.world.worldgen.cell.heightmap.ControlPoints;
import com.regenerationforrged.world.worldgen.noise.NoiseUtil;
import com.regenerationforrged.world.worldgen.noise.NoiseUtil.Vec2f;
import com.regenerationforrged.world.worldgen.noise.module.Line;

public class Island {
	private int id;
	private Segment[] segments;
	private float coast2;
	private float deepOcean2;
	private float shallowOcean;
	private float shallowOcean2;
	private float radius;
	private float deepMod;
	private float shallowMod;
	private Vec2f center;
	private Vec2f min;
	private Vec2f max;

	public Island(int id, Segment[] segments, ControlPoints controlPoints, float deepOcean, float shallowOcean, float coast, float inland) {
		float x = 0.0F;
		float y = 0.0F;
		int points = segments.length + 1;
		float minX = Float.MAX_VALUE;
		float minZ = Float.MAX_VALUE;
		float maxX = Float.MIN_VALUE;
		float maxZ = Float.MIN_VALUE;
		float maxRadius = coast;
		for (int i = 0; i < segments.length; ++i) {
			Segment segment = segments[i];
			minX = Math.min(minX, segment.minX());
			minZ = Math.min(minZ, segment.minY());
			maxX = Math.max(maxX, segment.maxX());
			maxZ = Math.max(maxZ, segment.maxY());
			maxRadius = Math.max(maxRadius, segment.maxScale() * coast);
			if (i == 0) {
				x += segment.a.x();
				y += segment.a.y();
			}
			x += segment.b.x();
			y += segment.b.y();
		}
		this.id = id;
		this.shallowOcean = shallowOcean;
		this.coast2 = coast * coast;
		this.deepOcean2 = deepOcean * deepOcean;
		this.shallowOcean2 = shallowOcean * shallowOcean;
		this.deepMod = 0.25F;
		this.shallowMod = 1.0F - this.deepMod;
		minX -= coast;
		minZ -= coast;
		maxX += coast;
		maxZ += coast;
		float maxDim = Math.max(maxX - minX, maxZ - minZ);
		this.center = new Vec2f(x / points, y / points);
		this.min = new Vec2f(minX - maxRadius, minZ - maxRadius);
		this.max = new Vec2f(maxX + maxRadius, maxZ + maxRadius);
		this.radius = maxDim;
	}

	public int getId() {
		return this.id;
	}

	public Segment[] getSegments() {
		return this.segments;
	}

	public float radius() {
		return this.radius;
	}

	public float coast() {
		return this.shallowOcean;
	}

	public void translate(Vec2f offset) {
		this.center = new Vec2f(this.center.x() + offset.x(), this.center.y() + offset.y());
		this.min = new Vec2f(this.min.x() + offset.x(), this.min.y() + offset.y());
		this.max = new Vec2f(this.max.x() + offset.x(), this.max.y() + offset.y());
		for (int i = 0; i < this.segments.length; ++i) {
			this.segments[i] = this.segments[i].translate(offset);
		}
	}

	public Vec2f getMin() {
		return this.min;
	}

	public Vec2f getMax() {
		return this.max;
	}

	public Vec2f getCenter() {
		return this.center;
	}

	public boolean overlaps(Island other) {
		return this.overlaps(other.min, other.max);
	}

	public boolean overlaps(Vec2f min, Vec2f max) {
		return this.min.x() < max.x() && this.max.x() > min.x() && this.min.y() < max.y() && this.max.y() > min.y();
	}

	public boolean contains(Vec2f vec) {
		return this.contains(vec.x(), vec.y());
	}

	public boolean contains(float x, float z) {
		return x > this.min.x() && x < this.max.x() && z > this.min.y() && z < this.max.y();
	}

	public float getEdgeValue(float x, float y) {
		float value = this.getEdgeDist2(x, y, this.deepOcean2);
		float deepValue = Math.min(this.deepOcean2, value);
		float shallowValue = Math.min(this.shallowOcean2, value);
		return this.process(deepValue, shallowValue);
	}

	public float getLandValue(float x, float y) {
		float value = this.getEdgeDist2(x, y, this.shallowOcean2);
		if (value < this.shallowOcean2) {
			value = (this.shallowOcean2 - value) / this.shallowOcean2;
			return NoiseUtil.curve(value, 0.75F, 4.0F);
		}
		return 0.0F;
	}

	private float getEdgeDist2(float x, float y, float minDist2) {
		float value = minDist2;
		for (Segment segment : this.segments) {
			float dx = segment.dx;
			float dy = segment.dy;
			float t = (x - segment.a.x()) * dx + (y - segment.a.y()) * dy;
			t /= segment.length2;
			float px;
			float py;
			float scale;
			if (t < 0.0F) {
				px = segment.a.x();
				py = segment.a.y();
				scale = segment.scale2A;
			} else if (t > 1.0F) {
				px = segment.b.x();
				py = segment.b.y();
				scale = segment.scale2B;
			} else {
				px = segment.a.x() + t * dx;
				py = segment.a.y() + t * dy;
				scale = NoiseUtil.lerp(segment.scale2A, segment.scale2B, t);
			}
			float v = Line.distSq(x, y, px, py) / scale;
			value = Math.min(v, value);
		}
		return value;
	}

	private float process(float deepValue, float shallowValue) {
		if (deepValue == this.deepOcean2) {
			return 0.0F;
		}
		if (deepValue > this.shallowOcean2) {
			deepValue = (deepValue - this.shallowOcean2) / (this.deepOcean2 - this.shallowOcean2);
			deepValue = 1.0F - deepValue;
			deepValue *= deepValue;
			return deepValue * this.deepMod;
		}
		if (shallowValue == this.shallowOcean2) {
			return this.deepMod;
		}
		if (shallowValue > this.coast2) {
			shallowValue = (shallowValue - this.coast2) / (this.shallowOcean2 - this.coast2);
			shallowValue = 1.0F - shallowValue;
			return this.deepMod + shallowValue * this.shallowMod;
		}
		return 1.0F;
	}
}

---

package com.regenerationforrged.world.worldgen.cell.continent.fancy;

import com.regenerationforrged.world.worldgen.noise.NoiseUtil.Vec2f;

public class Segment {
	public Vec2f a;
	public Vec2f b;
	public float dx;
	public float dy;
	public float length;
	public float length2;
	public float scaleA;
	public float scale2A;
	public float scaleB;
	public float scale2B;

	public Segment(Vec2f a, Vec2f b, float scaleA, float scaleB) {
		this.a = a;
		this.b = b;
		this.scaleA = scaleA;
		this.scaleB = scaleB;
		this.scale2A = scaleA * scaleA;
		this.scale2B = scaleB * scaleB;
		this.dx = b.x() - a.x();
		this.dy = b.y() - a.y();
		this.length = (float) Math.sqrt(this.dx * this.dx + this.dy * this.dy);
		this.length2 = this.dx * this.dx + this.dy * this.dy;
	}

	public float minX() {
		return Math.min(this.a.x(), this.b.x());
	}

	public float minY() {
		return Math.min(this.a.y(), this.b.y());
	}

	public float maxX() {
		return Math.max(this.a.x(), this.b.x());
	}

	public float maxY() {
		return Math.max(this.a.y(), this.b.y());
	}

	public float maxScale() {
		return Math.max(this.scaleA, this.scaleB);
	}

	public Segment translate(Vec2f offset) {
		return new Segment(new Vec2f(this.a.x() + offset.x(), this.a.y() + offset.y()), new Vec2f(this.b.x() + offset.x(), this.b.y() + offset.y()), this.scaleA, this.scaleB);
	}
}

---

package com.regenerationforrged.world.worldgen.cell.continent.infinite;

import com.regenerationforrged.world.worldgen.GeneratorContext;
import com.regenerationforrged.world.worldgen.cell.Cell;
import com.regenerationforrged.world.worldgen.cell.continent.SimpleContinent;
import com.regenerationforrged.world.worldgen.cell.continent.simple.SimpleRiverGenerator;
import com.regenerationforrged.world.worldgen.cell.rivermap.LegacyRiverCache;
import com.regenerationforrged.world.worldgen.cell.rivermap.RiverCache;
import com.regenerationforrged.world.worldgen.cell.rivermap.Rivermap;

public class InfiniteContinentGenerator implements SimpleContinent {
	private RiverCache riverCache;
	
	public InfiniteContinentGenerator(GeneratorContext context) {
        this.riverCache = new LegacyRiverCache(new SimpleRiverGenerator(this, context));
	}
	
	@Override
	public void apply(Cell cell, float x, float z) {
		cell.continentId = 0.0F;
		cell.continentEdge = 0.0F;
		cell.continentX = 0;
		cell.continentZ = 0;
	}

	@Override
	public Rivermap getRivermap(int x, int z) {
		return this.riverCache.getRivers(x, z);
	}

	@Override
	public long getNearestCenter(float x, float z) {
		return 0;
	}

	@Override
	public float getEdgeValue(float x, float z) {
		return 1.0F;
	}
}

---

package com.regenerationforrged.world.worldgen.cell.continent.simple;

import com.regenerationforrged.data.worldgen.preset.settings.WorldSettings;
import com.regenerationforrged.world.worldgen.GeneratorContext;
import com.regenerationforrged.world.worldgen.cell.Cell;
import com.regenerationforrged.world.worldgen.cell.continent.SimpleContinent;
import com.regenerationforrged.world.worldgen.cell.heightmap.ControlPoints;
import com.regenerationforrged.world.worldgen.cell.rivermap.LegacyRiverCache;
import com.regenerationforrged.world.worldgen.cell.rivermap.RiverCache;
import com.regenerationforrged.world.worldgen.cell.rivermap.Rivermap;
import com.regenerationforrged.world.worldgen.noise.NoiseUtil;
import com.regenerationforrged.world.worldgen.noise.NoiseUtil.Vec2f;
import com.regenerationforrged.world.worldgen.noise.domain.Domain;
import com.regenerationforrged.world.worldgen.noise.domain.Domains;
import com.regenerationforrged.world.worldgen.noise.function.DistanceFunction;
import com.regenerationforrged.world.worldgen.noise.function.EdgeFunction;
import com.regenerationforrged.world.worldgen.noise.module.Noise;
import com.regenerationforrged.world.worldgen.noise.module.Noises;
import com.regenerationforrged.world.worldgen.util.PosUtil;
import com.regenerationforrged.world.worldgen.util.Seed;

public abstract class ContinentGenerator implements SimpleContinent {
    protected int seed;
    protected float frequency;
    protected int continentScale;
    private DistanceFunction distanceFunc;
    private ControlPoints controlPoints;
    private float clampMin;
    private float clampMax;
    private float clampRange;
    private float offsetAlpha;
    protected Domain warp;
    protected Noise shape;
    protected RiverCache cache;
    
    public ContinentGenerator(Seed seed, GeneratorContext context) {
        WorldSettings settings = context.preset.world();
        int tectonicScale = settings.continent.continentScale * 4;
        this.continentScale = settings.continent.continentScale / 2;
        this.seed = seed.next();
        this.distanceFunc = settings.continent.continentShape;
        this.controlPoints = ControlPoints.make(settings.controlPoints);
        this.frequency = 1.0F / tectonicScale;
        this.clampMin = 0.2F;
        this.clampMax = 1.0F;
        this.clampRange = this.clampMax - this.clampMin;
        this.offsetAlpha = context.preset.world().continent.continentJitter;
        
        Domain warp = Domains.domainPerlin(seed.next(), 20, 2, 20.0F);
        warp = Domains.compound(warp, Domains.domainSimplex(seed.next(), this.continentScale, 3, this.continentScale));
        this.warp = warp;
        
        Noise shape = Noises.simplex(seed.next(), settings.continent.continentScale * 2, 1);
        shape = Noises.add(shape, 0.65F);
        shape = Noises.clamp(shape, 0.0F, 1.0F);
        this.shape = shape;

        this.cache = new LegacyRiverCache(new SimpleRiverGenerator(this, context));
    }
    
    @Override
    public Rivermap getRivermap(int x, int y) {
        return this.cache.getRivers(x, y);
    }
    
    @Override
    public void apply(Cell cell, float x, float y) {
        float ox = this.warp.getOffsetX(x, y, 0);
        float oz = this.warp.getOffsetZ(x, y, 0);
        float px = x + ox;
        float py = y + oz;
        px *= this.frequency;
        py *= this.frequency;
        int xr = NoiseUtil.floor(px);
        int yr = NoiseUtil.floor(py);
        int cellX = xr;
        int cellY = yr;
        float centerX = px;
        float centerY = py;
        float edgeDistance = 999999.0F;
        float edgeDistance2 = 999999.0F;
        for (int dy = -1; dy <= 1; ++dy) {
            for (int dx = -1; dx <= 1; ++dx) {
                int cx = xr + dx;
                int cy = yr + dy;
                Vec2f vec = NoiseUtil.cell(this.seed, cx, cy);
                float cxf = cx + vec.x() * this.offsetAlpha;
                float cyf = cy + vec.y() * this.offsetAlpha;
                float distance = this.distanceFunc.apply(cxf - px, cyf - py);
                if (distance < edgeDistance) {
                    edgeDistance2 = edgeDistance;
                    edgeDistance = distance;
                    centerX = cxf;
                    centerY = cyf;
                    cellX = cx;
                    cellY = cy;
                }
                else if (distance < edgeDistance2) {
                    edgeDistance2 = distance;
                }
            }
        }
        cell.continentId = this.cellIdentity(this.seed, cellX, cellY);
        cell.continentEdge = this.cellEdgeValue(edgeDistance, edgeDistance2);
        cell.continentX = (int) (centerX / this.frequency);
        cell.continentZ = (int) (centerY / this.frequency);
        cell.continentEdge *= this.getShape(x, y, cell.continentEdge);
    }
    
    @Override
    public float getEdgeValue(float x, float y) {
        float ox = this.warp.getOffsetX(x, y, 0);
        float oz = this.warp.getOffsetZ(x, y, 0);
        float px = x + ox;
        float py = y + oz;
        px *= this.frequency;
        py *= this.frequency;
        int xr = NoiseUtil.floor(px);
        int yr = NoiseUtil.floor(py);
        float edgeDistance = 999999.0F;
        float edgeDistance2 = 999999.0F;
        for (int dy = -1; dy <= 1; ++dy) {
            for (int dx = -1; dx <= 1; ++dx) {
                int cx = xr + dx;
                int cy = yr + dy;
                Vec2f vec = NoiseUtil.cell(this.seed, cx, cy);
                float cxf = cx + vec.x() * this.offsetAlpha;
                float cyf = cy + vec.y() * this.offsetAlpha;
                float distance = this.distanceFunc.apply(cxf - px, cyf - py);
                if (distance < edgeDistance) {
                    edgeDistance2 = edgeDistance;
                    edgeDistance = distance;
                }
                else if (distance < edgeDistance2) {
                    edgeDistance2 = distance;
                }
            }
        }
        float edgeValue = this.cellEdgeValue(edgeDistance, edgeDistance2);
        float shapeNoise = this.getShape(x, y, edgeValue);
        return edgeValue * shapeNoise;
    }
    
    @Override
    public long getNearestCenter(float x, float z) {
        float ox = this.warp.getOffsetX(x, z, 0);
        float oz = this.warp.getOffsetZ(x, z, 0);
        float px = x + ox;
        float py = z + oz;
        px *= this.frequency;
        py *= this.frequency;
        float centerX = px;
        float centerY = py;
        int xr = NoiseUtil.floor(px);
        int yr = NoiseUtil.floor(py);
        float edgeDistance = 999999.0F;
        for (int dy = -1; dy <= 1; ++dy) {
            for (int dx = -1; dx <= 1; ++dx) {
                int cx = xr + dx;
                int cy = yr + dy;
                Vec2f vec = NoiseUtil.cell(this.seed, cx, cy);
                float cxf = cx + vec.x() * this.offsetAlpha;
                float cyf = cy + vec.y() * this.offsetAlpha;
                float distance = this.distanceFunc.apply(cxf - px, cyf - py);
                if (distance < edgeDistance) {
                    edgeDistance = distance;
                    centerX = cxf;
                    centerY = cyf;
                }
            }
        }
        int conX = (int)(centerX / this.frequency);
        int conZ = (int)(centerY / this.frequency);
        return PosUtil.pack(conX, conZ);
    }
    
    @Override
    public float getDistanceToOcean(int cx, int cz, float dx, float dz) {
        float high = this.getDistanceToEdge(cx, cz, dx, dz);
        float low = 0.0F;
        for (int i = 0; i < 50; ++i) {
            float mid = (low + high) / 2.0F;
            float x = cx + dx * mid;
            float z = cz + dz * mid;
            float edge = this.getEdgeValue(x, z);
            if (edge > this.controlPoints.shallowOcean()) {
                low = mid;
            }
            else {
                high = mid;
            }
            if (high - low < 10.0f) {
                break;
            }
        }
        return high;
    }
    
    @Override
    public float getDistanceToEdge(int cx, int cz, float dx, float dz) {
        float distance = (float)(this.continentScale * 4);
        for (int i = 0; i < 10; ++i) {
            float x = cx + dx * distance;
            float z = cz + dz * distance;
            long centerPos = this.getNearestCenter(x, z);
            int conX = PosUtil.unpackLeft(centerPos);
            int conZ = PosUtil.unpackRight(centerPos);
            distance += distance;
            if (conX != cx || conZ != cz) {
                float low = 0.0F;
                float high = distance;
                for (int j = 0; j < 50; ++j) {
                    float mid = (low + high) / 2.0F;
                    float px = cx + dx * mid;
                    float pz = cz + dz * mid;
                    centerPos = this.getNearestCenter(px, pz);
                    conX = PosUtil.unpackLeft(centerPos);
                    conZ = PosUtil.unpackRight(centerPos);
                    if (conX == cx && conZ == cz) {
                        low = mid;
                    }
                    else {
                        high = mid;
                    }
                    if (high - low < 50.0F) {
                        break;
                    }
                }
                return high;
            }
        }
        return distance;
    }
    
    protected float cellIdentity(int seed, int cellX, int cellY) {
        float value = NoiseUtil.valCoord2D(seed, cellX, cellY);
        return NoiseUtil.map(value, -1.0F, 1.0F, 2.0F);
    }
    
    protected float cellEdgeValue(float distance, float distance2) {
    	EdgeFunction edge = EdgeFunction.DISTANCE_2_DIV;
        float value = edge.apply(distance, distance2);
        value = 1.0F - NoiseUtil.map(value, edge.min(), edge.max(), edge.range());
        if (value <= this.clampMin) {
            return 0.0F;
        }
        if (value >= this.clampMax) {
            return 1.0F;
        }
        return (value - this.clampMin) / this.clampRange;
    }
    
    protected float getShape(float x, float z, float edgeValue) {
        if (edgeValue >= this.controlPoints.inland()) {
            return 1.0F;
        }
        float alpha = edgeValue / this.controlPoints.inland();
        return this.shape.compute(x, z, 0) * alpha;
    }
}

---


package com.regenerationforrged.world.worldgen.cell.continent.simple;

import com.regenerationforrged.world.worldgen.GeneratorContext;
import com.regenerationforrged.world.worldgen.util.Seed;

public class MultiContinentGenerator extends ContinentGenerator {
	
    public MultiContinentGenerator(Seed seed, GeneratorContext context) {
        super(seed, context);
    }
}

---


package com.regenerationforrged.world.worldgen.cell.continent.simple;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import com.regenerationforrged.world.worldgen.GeneratorContext;
import com.regenerationforrged.world.worldgen.cell.continent.SimpleContinent;
import com.regenerationforrged.world.worldgen.cell.rivermap.gen.GenWarp;
import com.regenerationforrged.world.worldgen.cell.rivermap.river.BaseRiverGenerator;
import com.regenerationforrged.world.worldgen.cell.rivermap.river.Network;
import com.regenerationforrged.world.worldgen.cell.rivermap.river.River;
import com.regenerationforrged.world.worldgen.cell.rivermap.river.RiverCarver;
import com.regenerationforrged.world.worldgen.cell.rivermap.river.RiverWarp;
import com.regenerationforrged.world.worldgen.noise.NoiseUtil;

public class SimpleRiverGenerator extends BaseRiverGenerator<SimpleContinent> {
	
	public SimpleRiverGenerator(SimpleContinent continent, GeneratorContext context) {
		super(continent, context);
	}

	@Override
	public List<Network.Builder> generateRoots(int x, int z, Random random, GenWarp warp) {
		float start = random.nextFloat();
		float spacing = 6.2831855F / this.count;
		float spaceVar = spacing * 0.75F;
		float spaceBias = -spaceVar / 2.0F;
		List<Network.Builder> roots = new ArrayList<>(this.count);
		for (int i = 0; i < this.count; ++i) {
			float variance = random.nextFloat() * spaceVar + spaceBias;
			float angle = start + spacing * i + variance;
			float dx = NoiseUtil.sin(angle);
			float dz = NoiseUtil.cos(angle);
			float startMod = 0.05F + random.nextFloat() * 0.45F;
			float length = this.continent.getDistanceToOcean(x, z, dx, dz);
			float startDist = Math.max(400.0F, startMod * length);
			float x2 = x + dx * startDist;
			float z2 = z + dz * startDist;
			float x3 = x + dx * length;
			float z3 = z + dz * length;
			float valleyWidth = 275.0F * River.MAIN_VALLEY.next(random);
			River river = new River((float) (int) x2, (float) (int) z2, (float) (int) x3, (float) (int) z3);
			RiverCarver.Settings settings = BaseRiverGenerator.creatSettings(random);
			settings.fadeIn = this.main.fade;
			settings.valleySize = valleyWidth;
			RiverWarp riverWarp = RiverWarp.create(0.1F, 0.85F, random);
			RiverCarver carver = new RiverCarver(river, riverWarp, this.main, settings, this.levels);
			Network.Builder branch = Network.builder(carver);
			roots.add(branch);
			this.addLake(branch, random, warp);
		}
		return roots;
	}
}

---


package com.regenerationforrged.world.worldgen.cell.continent.simple;

import com.regenerationforrged.world.worldgen.GeneratorContext;
import com.regenerationforrged.world.worldgen.cell.Cell;
import com.regenerationforrged.world.worldgen.noise.NoiseUtil.Vec2i;
import com.regenerationforrged.world.worldgen.util.PosUtil;
import com.regenerationforrged.world.worldgen.util.Seed;

public class SingleContinentGenerator extends ContinentGenerator {
    private Vec2i center;
    
    public SingleContinentGenerator(Seed seed, GeneratorContext context) {
        super(seed, context);
        long center = this.getNearestCenter(0.0F, 0.0F);
        int cx = PosUtil.unpackLeft(center);
        int cz = PosUtil.unpackRight(center);
        this.center = new Vec2i(cx, cz);
    }
    
    @Override
    public void apply(Cell cell, float x, float y) {
        super.apply(cell, x, y);
        if (cell.continentX != this.center.x() || cell.continentZ != this.center.y()) {
            cell.continentId = 0.0F;
            cell.continentEdge = 0.0F;
            cell.continentX = 0;
            cell.continentZ = 0;
        }
    }
}
