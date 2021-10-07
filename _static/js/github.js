const githubOrganization = 'FreeFem'
const githubRepository = 'FreeFem-sources'

function HTTPGet(url, callback) {
   const xhr = new XMLHttpRequest()
   xhr.open('GET', url, true)
   xhr.onload = function(e) {
      if (xhr.readyState === 4) {
         if (xhr.status === 200) {
            callback(xhr.responseText)
         } else {
            console.error(xhr.statusText)
         }
      }
   }
   xhr.onerror = function(e) {
      console.error(xhr.statusText)
   }
   xhr.send(null)
}

const githubURL = 'https://api.github.com/repos/' + githubOrganization + '/' + githubRepository

function starsAndForks(data) {
   const jdata = JSON.parse(data)

   const headerGithubStarsForks = document.getElementById('headerGithubStarsForks')
   headerGithubStarsForks.style.display = 'block'

   const stars = jdata.stargazers_count
   const headerGithubStars = document.getElementById('headerGithubStars')
   headerGithubStars.innerHTML = stars

   const forks = jdata.forks_count
   const headerGithubForks = document.getElementById('headerGithubForks')
   headerGithubForks.innerHTML = forks
}

HTTPGet(githubURL, starsAndForks)
