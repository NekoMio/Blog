import { getAssetFromKV, mapRequestToAsset } from '@cloudflare/kv-asset-handler'

/**
 * The DEBUG flag will do two things that help during development:
 * 1. we will skip caching on the edge, which makes it easier to
 *    debug.
 * 2. we will return an error message on exception in your Response rather
 *    than the default 404.html page.
 */
const DEBUG = false

addEventListener('fetch', event => {
  try {
    event.respondWith(handleEvent(event))
  } catch (e) {
    if (DEBUG) {
      return event.respondWith(
        new Response(e.message || e.toString(), {
          status: 500,
        }),
      )
    }
    event.respondWith(new Response('Internal Error', { status: 500 }))
  }
})

async function handleEvent(event) {
  const { origin, pathname: path, search } = new URL(event.request.url);
  if (path.endsWith('/index.html')) {
    return new Response(null, {
      status: 301,
        headers: {
          'Location': `${origin}${path.substring(0, path.length - 10)}${search}`,
          'Cache-Control': 'max-age=3600'
        }
    });
  }
  if (path === '/atom.xml') {
    return getAssetFromKV(event, {
      cacheControl: {
        edgeTtl: 60 * 60,
        browserTtl: 2 * 60 * 60
      }
    });
  }
  
  if (path.startsWith('/css/'))  {
    const response = await getAssetFromKV(event, {
      cacheControl: {
        edgeTtl: 365 * 24 * 60 * 60,
        browserTtl: 365 * 24 * 60 * 60
      }
    });
    // getAssetFromKV 不会添加 immutable，所以需要手动覆盖掉 Cache-Control
    response.headers.set('cache-control', `public, max-age=${365 * 24 * 60 * 60}, immutable`);
    return response;
  }
  
  const response = await getAssetFromKV(event, {
    cacheControl: {
      edgeTtl: 60 * 60,
      browserTtl: 5 * 60
    }
  });
  
  response.headers.set('X-XSS-Protection', '1; mode=block');
  return response;
}